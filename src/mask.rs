//! Q2 masking of M-biased bases.
//!
//! Once mask lengths are frozen from the M-bias curves, each primary mapped
//! record has its biased qualities set to a low value (`--mbias-mask-quality`,
//! default 2) so downstream base-quality-aware callers ignore them. Nothing
//! else about the record changes — no clip, no POS/CIGAR/tag/mate rewrite — and
//! masked bases fall below methylsieve's own `--min-base-quality` gate, so they
//! drop out of its tallies for free.
//!
//! What gets masked, per primary mapped record:
//! - **Own 5':** the first `K` sequencing cycles (low stored positions for a
//!   forward read, high for a reverse read — SEQ is stored genomic-forward).
//! - **Single-end 3':** additionally the last `K_3'` cycles (its far template
//!   end is unknown, so both ends are learned and masked).
//! - **Orphan 3' mirror:** a paired read whose mate is unmapped masks its 3' end
//!   by the *mate role's* 5' length, in case it read through the whole template
//!   (no mate is present to propagate from).
//! - **Shared masked positions (any pair orientation):** a reference position
//!   masked in one mate is also masked in the other wherever that mate covers it.
//!   This is the masking analogue of the overlap "count each position once" rule,
//!   driven by reference coverage rather than strand — so it handles FR, FF, RF,
//!   and dovetailed pairs alike, and never touches a mate that doesn't actually
//!   cover the masked position.

use fgumi_raw_bam::RawRecord;

use crate::mbias::{DetectParams, MbiasAccumulator, ReadEnd, ReadRole, detect_mask_length};
use crate::record::{
    FLAG_FIRST_SEGMENT, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_REVERSE, has, is_primary_mapped,
    read_role, ref_span_for_query_window,
};

/// Frozen 5'/3' mask lengths (in sequencing cycles) per read role, plus the
/// quality value masked positions are set to.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MaskPlan {
    k_r1_5p: usize,
    k_r2_5p: usize,
    k_se_5p: usize,
    k_se_3p: usize,
    mask_quality: u8,
}

impl MaskPlan {
    /// Learn the plan from accumulated M-bias (CpG curves).
    #[must_use]
    pub(crate) fn learn(acc: &MbiasAccumulator, detect: DetectParams, mask_quality: u8) -> Self {
        let d = |role, end| detect_mask_length(acc, role, end, detect);
        Self {
            k_r1_5p: d(ReadRole::R1, ReadEnd::FivePrime),
            k_r2_5p: d(ReadRole::R2, ReadEnd::FivePrime),
            k_se_5p: d(ReadRole::Se, ReadEnd::FivePrime),
            k_se_3p: d(ReadRole::Se, ReadEnd::ThreePrime),
            mask_quality,
        }
    }

    /// Construct an explicit plan with given R1/R2 5' lengths (tests only).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn explicit(k_r1_5p: usize, k_r2_5p: usize, mask_quality: u8) -> Self {
        Self { k_r1_5p, k_r2_5p, k_se_5p: 0, k_se_3p: 0, mask_quality }
    }

    /// The 5' mask length for a role.
    fn role_5p(&self, role: ReadRole) -> usize {
        match role {
            ReadRole::R1 => self.k_r1_5p,
            ReadRole::R2 => self.k_r2_5p,
            ReadRole::Se => self.k_se_5p,
        }
    }

    /// One-line summary of the frozen lengths, for the run log.
    #[must_use]
    pub(crate) fn summary(&self) -> String {
        format!(
            "R1 5'={}, R2 5'={}, SE 5'={}, SE 3'={} (mask Q{})",
            self.k_r1_5p, self.k_r2_5p, self.k_se_5p, self.k_se_3p, self.mask_quality
        )
    }
}

/// Apply `plan` to every primary mapped record of one template, in place.
///
/// Three phases so cross-mate propagation can read both mates before any
/// mutation: (1) each record's *own* mask windows (5', plus SE 3' or orphan
/// mirror); (2) propagate masked reference positions between the two primary
/// mates by coverage; (3) write Q2 over all collected windows.
pub(crate) fn mask_template(plan: &MaskPlan, recs: &mut [RawRecord]) {
    let n = recs.len();

    // Locate the two primary mapped paired mates (for cross-mate propagation).
    let (mut r1, mut r2) = (None, None);
    for (i, rec) in recs.iter().enumerate() {
        let f = rec.flags();
        if !is_primary_mapped(f) || !has(f, FLAG_PAIRED) {
            continue;
        }
        if has(f, FLAG_FIRST_SEGMENT) {
            r1.get_or_insert(i);
        } else {
            r2.get_or_insert(i);
        }
    }

    // Phase 1 — each primary mapped record's own mask windows (stored positions).
    let mut windows: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for i in 0..n {
        let f = recs[i].flags();
        if !is_primary_mapped(f) {
            continue;
        }
        let seq_len = recs[i].l_seq() as usize;
        if seq_len == 0 {
            continue;
        }
        let reverse = has(f, FLAG_REVERSE);
        let role = read_role(f);

        // Own 5': low stored for a forward read, high for a reverse read.
        let k5 = plan.role_5p(role).min(seq_len);
        if k5 > 0 {
            windows[i].push(if reverse { (seq_len - k5, seq_len) } else { (0, k5) });
        }

        if role == ReadRole::Se {
            // Single-end: also the learned 3' end.
            let k3 = plan.k_se_3p.min(seq_len);
            if k3 > 0 {
                windows[i].push(if reverse { (0, k3) } else { (seq_len - k3, seq_len) });
            }
        } else {
            // Orphan (mate unmapped or absent): no mate to propagate from, so
            // mirror the mate role's 5' onto the 3' end in case the read ran
            // through the whole template. Present mates are handled in phase 2.
            let mate = if has(f, FLAG_FIRST_SEGMENT) { r2 } else { r1 };
            if has(f, FLAG_MATE_UNMAPPED) || mate.is_none() {
                let km = plan.role_5p(role.mate()).min(seq_len);
                if km > 0 {
                    windows[i].push(if reverse { (0, km) } else { (seq_len - km, seq_len) });
                }
            }
        }
    }

    // Phase 2 — propagate masked reference positions between the two primary
    // mates (any orientation), when both are present, mapped, and co-located.
    if let (Some(i1), Some(i2)) = (r1, r2)
        && recs[i1].ref_id() >= 0
        && recs[i1].ref_id() == recs[i2].ref_id()
    {
        let masked_ref_ranges = |idx: usize| -> Vec<(usize, usize)> {
            windows[idx]
                .iter()
                .filter_map(|&(lo, hi)| ref_span_for_query_window(&recs[idx], lo, hi))
                .collect()
        };
        let r1_ranges = masked_ref_ranges(i1);
        let r2_ranges = masked_ref_ranges(i2);
        // Whatever R2 masked, mask in R1 where R1 covers it, and vice versa.
        for &(lo, hi) in &r2_ranges {
            if let Some(w) = stored_window_for_ref(&recs[i1], lo, hi) {
                windows[i1].push(w);
            }
        }
        for &(lo, hi) in &r1_ranges {
            if let Some(w) = stored_window_for_ref(&recs[i2], lo, hi) {
                windows[i2].push(w);
            }
        }
    }

    // Phase 3 — write Q2 over every collected window.
    for i in 0..n {
        if windows[i].is_empty() {
            continue;
        }
        let q = recs[i].quality_scores_mut();
        for &(lo, hi) in &windows[i] {
            let hi = hi.min(q.len());
            for b in &mut q[lo.min(hi)..hi] {
                *b = plan.mask_quality;
            }
        }
    }
}

/// Stored half-open window `[q_lo, q_hi)` whose aligned reference positions fall
/// in `[ref_lo, ref_hi)`. The reference-to-query inverse of
/// [`ref_span_for_query_window`]. `None` when no aligned base maps into that
/// reference range.
fn stored_window_for_ref(rec: &RawRecord, ref_lo: usize, ref_hi: usize) -> Option<(usize, usize)> {
    if ref_lo >= ref_hi {
        return None;
    }
    let mut qpos = 0usize;
    let mut rpos = rec.pos().max(0) as usize;
    let (mut lo, mut hi) = (usize::MAX, 0usize);
    for op in rec.cigar_ops_iter() {
        let len = (op >> 4) as usize;
        match op & 0xf {
            0 | 7 | 8 => {
                let a = rpos.max(ref_lo);
                let b = (rpos + len).min(ref_hi);
                if a < b {
                    lo = lo.min(qpos + (a - rpos));
                    hi = hi.max(qpos + (b - rpos));
                }
                qpos += len;
                rpos += len;
            }
            1 | 4 => qpos += len,
            2 | 3 => rpos += len,
            _ => {}
        }
        if rpos >= ref_hi {
            break;
        }
    }
    (lo < hi).then_some((lo, hi))
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, Cursor};

    use super::*;
    use crate::record::FLAG_LAST_SEGMENT;
    use crate::sam_reader::SamReader;

    /// Parse SAM lines into `RawRecord`s (correct BAM layout).
    fn parse(records: &[&str], contig_len: usize) -> Vec<RawRecord> {
        let mut sam = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{contig_len}\n");
        for r in records {
            sam.push_str(r);
            sam.push('\n');
        }
        let boxed: Box<dyn BufRead> = Box::new(Cursor::new(sam.into_bytes()));
        let mut reader = SamReader::new(boxed);
        reader.read_header().unwrap();
        let mut out = Vec::new();
        let mut rec = RawRecord::new();
        while reader.read_record(&mut rec).unwrap() {
            out.push(std::mem::replace(&mut rec, RawRecord::new()));
        }
        out
    }

    fn line(qname: &str, flag: u16, pos: u32, cigar: &str, seq: &str) -> String {
        let qual = "I".repeat(seq.len());
        format!("{qname}\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
    }

    /// Quality bytes of a record (Phred, not ASCII). Phred 40 = input 'I'.
    fn quals(rec: &RawRecord) -> Vec<u8> {
        rec.quality_scores().to_vec()
    }

    /// Flags for a lone primary R1 (paired bit set, no mate in template → the
    /// orphan path, with k_r2 = 0 here so only the 5' mask applies).
    const R1: u16 = FLAG_PAIRED | FLAG_FIRST_SEGMENT;

    #[test]
    fn forward_read_masks_low_stored_5p() {
        let mut recs = parse(&[&line("r", R1, 1, "10M", "CACACACACA")], 20);
        mask_template(&MaskPlan::explicit(3, 0, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[0..3], &[2, 2, 2], "first 3 cycles masked");
        assert!(q[3..].iter().all(|&b| b == 40), "rest untouched");
    }

    #[test]
    fn reverse_read_masks_high_stored_5p() {
        // Reverse read: 5' end is the HIGH stored position.
        let mut recs = parse(&[&line("r", R1 | FLAG_REVERSE, 1, "10M", "CACACACACA")], 20);
        mask_template(&MaskPlan::explicit(3, 0, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[7..10], &[2, 2, 2], "last 3 cycles masked (5' of reverse read)");
        assert!(q[..7].iter().all(|&b| b == 40));
    }

    #[test]
    fn read_shorter_than_k_is_fully_masked() {
        let mut recs = parse(&[&line("r", R1, 1, "4M", "CACA")], 20);
        mask_template(&MaskPlan::explicit(10, 0, 2), &mut recs);
        assert!(quals(&recs[0]).iter().all(|&b| b == 2), "all masked when K > read length");
    }

    /// FR pair, fully overlapping. R1 masks its 5' (ref [0,3)); R2 (reverse)
    /// masks its 5' (ref [6,10)). Propagation masks each mate at the OTHER's
    /// masked reference positions, so both end up masked at [0,3) and [6,10).
    #[test]
    fn mate_propagation_masks_shared_reference_positions() {
        let r1 = line("p", FLAG_PAIRED | FLAG_FIRST_SEGMENT, 1, "10M", "CACACACACA");
        let r2 = line("p", FLAG_PAIRED | FLAG_LAST_SEGMENT | FLAG_REVERSE, 1, "10M", "CACACACACA");
        let mut recs = parse(&[&r1, &r2], 20);
        mask_template(&MaskPlan::explicit(3, 4, 2), &mut recs); // k_r1=3, k_r2=4
        let q1 = quals(&recs[0]); // R1 forward: own [0,3), propagated [6,10)
        assert_eq!(&q1[0..3], &[2, 2, 2], "R1 own 5'");
        assert!(q1[3..6].iter().all(|&b| b == 40), "R1 middle kept");
        assert_eq!(&q1[6..10], &[2, 2, 2, 2], "R1 masked at R2's masked positions");
        let q2 = quals(&recs[1]); // R2 reverse: own [6,10), propagated [0,3)
        assert_eq!(&q2[0..3], &[2, 2, 2], "R2 masked at R1's masked positions");
        assert!(q2[3..6].iter().all(|&b| b == 40), "R2 middle kept");
        assert_eq!(&q2[6..10], &[2, 2, 2, 2], "R2 own 5'");
    }

    /// Non-overlapping FR pair: the mates share no reference positions, so NO
    /// propagation occurs — only each read's own 5' is masked. (This is the case
    /// the old strand-proxy overhang logic wrongly masked in full.)
    #[test]
    fn non_overlapping_pair_does_not_propagate() {
        let r1 = line("p", FLAG_PAIRED | FLAG_FIRST_SEGMENT, 1, "10M", "CACACACACA");
        // R2 far downstream at ref [20,30); shares nothing with R1's [0,10).
        let r2 = line("p", FLAG_PAIRED | FLAG_LAST_SEGMENT | FLAG_REVERSE, 21, "10M", "CACACACACA");
        let mut recs = parse(&[&r1, &r2], 40);
        mask_template(&MaskPlan::explicit(3, 4, 2), &mut recs);
        let q1 = quals(&recs[0]);
        assert_eq!(&q1[0..3], &[2, 2, 2], "R1 own 5' only");
        assert!(q1[3..].iter().all(|&b| b == 40), "no spurious propagation into R1");
        let q2 = quals(&recs[1]);
        assert_eq!(&q2[6..10], &[2, 2, 2, 2], "R2 own 5' only");
        assert!(q2[..6].iter().all(|&b| b == 40), "no spurious propagation into R2");
    }

    /// Propagation is orientation-agnostic: a same-strand (FF) overlapping pair
    /// still propagates by reference coverage. R1 masks ref [0,3); R2 (k=0, no own
    /// mask) gets [0,3) masked purely from propagation.
    #[test]
    fn same_strand_pair_propagates_by_coverage() {
        let r1 = line("p", FLAG_PAIRED | FLAG_FIRST_SEGMENT, 1, "10M", "CACACACACA");
        let r2 = line("p", FLAG_PAIRED | FLAG_LAST_SEGMENT, 1, "10M", "CACACACACA"); // forward
        let mut recs = parse(&[&r1, &r2], 20);
        mask_template(&MaskPlan::explicit(3, 0, 2), &mut recs); // k_r2 = 0
        let q2 = quals(&recs[1]);
        assert_eq!(&q2[0..3], &[2, 2, 2], "R2 masked at R1's masked positions despite k_r2=0");
        assert!(q2[3..].iter().all(|&b| b == 40));
    }

    #[test]
    fn orphan_mirrors_mate_length_on_3p() {
        // R1 forward, mate unmapped (orphan). Own 5' = k_r1 (2); 3' mirror = k_r2 (3).
        let r1 = line(
            "o",
            FLAG_PAIRED | FLAG_FIRST_SEGMENT | FLAG_MATE_UNMAPPED,
            1,
            "10M",
            "CACACACACA",
        );
        let mut recs = parse(&[&r1], 20);
        mask_template(&MaskPlan::explicit(2, 3, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[0..2], &[2, 2], "own 5' (k_r1) masked");
        assert_eq!(&q[7..10], &[2, 2, 2], "3' mirrored by mate role k_r2");
        assert!(q[2..7].iter().all(|&b| b == 40));
    }
}
