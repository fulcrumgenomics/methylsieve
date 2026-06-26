//! Q2 masking of M-biased bases.
//!
//! Once mask lengths are frozen from the M-bias curves, each maskable record has
//! its biased qualities set to a low value (`--mbias-mask-quality`, default 2) so
//! downstream base-quality-aware callers ignore them. Nothing else about the
//! record changes — no clip, no POS/CIGAR/tag/mate rewrite.
//!
//! The mask is computed as a set of stored-position windows
//! ([`compute_mask_windows`]) and applied two ways from that one geometry: the Q2
//! write ([`apply_windows`]) for downstream callers, and — when masking would
//! drop a base below `--min-base-quality` — an exclusion from methylsieve's own
//! per-template tally and unconverted decision (see `RecordSkips` in
//! [`crate::sieve`]). Sharing one geometry keeps the reported metrics and the
//! call describing exactly the bases the masked output presents.
//!
//! **Which records are masked.** Every record that has SEQ and real base
//! qualities and is **not** a secondary alignment — i.e. primary, supplementary,
//! and even unmapped records. Secondary alignments are *never* masked: their SEQ
//! is frequently absent or a hard-clipped duplicate of the primary, and
//! base-quality-aware methylation callers ignore secondaries, so masking them is
//! pure risk for no downstream benefit. (M-bias *measurement* — the mask lengths
//! themselves — is learned from primary mapped records only, **before** masking,
//! so the curve is a pre-mask view of the data; only the tally and decision act
//! on the masked result.)
//!
//! What gets masked, per maskable record:
//! - **Own 5':** the first `K` sequencing cycles (low stored positions for a
//!   forward read, high for a reverse read — SEQ is stored genomic-forward).
//!   Cycles hard-clipped off the 5' end (e.g. a split read's supplementary) are
//!   absent from SEQ, so the stored window is shifted by the hard-clip length and
//!   may be empty.
//! - **Single-end 3':** additionally the last `K_3'` cycles (its far template
//!   end is unknown, so both ends are learned and masked).
//! - **Mate-5' mirror onto 3':** a paired read that can't recover its mate's
//!   masked positions from reference coverage — because the read itself is
//!   unmapped, or its mate is unmapped/absent — masks its 3' end by the *mate
//!   role's* 5' length, in case it read through to the mate's 5'. (A defensive
//!   over-estimate; harmless, since such reads carry no methylation calls unless
//!   later realigned.)
//! - **Shared masked positions (any orientation, any contig):** a reference
//!   position masked in one mate is masked in the other wherever that mate covers
//!   it, matched by contig. This is the masking analogue of the overlap "count
//!   each position once" rule, driven by reference coverage rather than strand —
//!   so it handles FR/FF/RF/dovetailed pairs and split-read supplementaries that
//!   cross an SV breakpoint onto the mate's contig alike, and never touches a
//!   mate that doesn't actually cover the masked position.

use fgumi_raw_bam::RawRecord;
use smallvec::SmallVec;

use crate::mbias::{
    DetectParams, MbiasAccumulator, ReadEnd, ReadRole, cpg_plateau, detect_mask_length,
};
use crate::record::{
    FLAG_FIRST_SEGMENT, FLAG_LAST_SEGMENT, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_REVERSE,
    FLAG_SECONDARY, FLAG_UNMAPPED, has, read_role, ref_span_for_query_window,
};

/// Stored-position mask windows for a single record. Inline-allocated: a record
/// almost always has at most its own 5' window plus one propagated overhang, so
/// the common case never touches the heap.
pub(crate) type MaskWindows = SmallVec<[(usize, usize); 2]>;

/// Frozen 5'/3' mask lengths (in sequencing cycles) per read role, plus the
/// quality value masked positions are set to.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MaskPlan {
    k_r1_5p: usize,
    k_r2_5p: usize,
    k_se_5p: usize,
    k_se_3p: usize,
    mask_quality: u8,
    /// CpG methylation plateau per role (by [`ReadRole::index`]) for the 5' end,
    /// `None` where no plateau was estimable. Combined with `plateau_fraction`
    /// this is the detection threshold the masks were cut at — exposed for plots.
    plateau_5p: [Option<f64>; 3],
    /// The `--mbias-plateau-fraction` the plan was learned with.
    plateau_fraction: f64,
}

impl MaskPlan {
    /// Learn the plan from accumulated M-bias (CpG curves).
    #[must_use]
    pub(crate) fn learn(acc: &MbiasAccumulator, detect: DetectParams, mask_quality: u8) -> Self {
        let d = |role, end| detect_mask_length(acc, role, end, detect);
        let mut plateau_5p = [None; 3];
        for role in ReadRole::ALL {
            plateau_5p[role.index()] = cpg_plateau(acc, role, ReadEnd::FivePrime, detect);
        }
        Self {
            k_r1_5p: d(ReadRole::R1, ReadEnd::FivePrime),
            k_r2_5p: d(ReadRole::R2, ReadEnd::FivePrime),
            k_se_5p: d(ReadRole::Se, ReadEnd::FivePrime),
            k_se_3p: d(ReadRole::Se, ReadEnd::ThreePrime),
            mask_quality,
            plateau_5p,
            plateau_fraction: detect.plateau_fraction,
        }
    }

    /// Construct an explicit plan with given R1/R2 5' lengths (tests only).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn explicit(k_r1_5p: usize, k_r2_5p: usize, mask_quality: u8) -> Self {
        Self {
            k_r1_5p,
            k_r2_5p,
            k_se_5p: 0,
            k_se_3p: 0,
            mask_quality,
            plateau_5p: [None; 3],
            plateau_fraction: DetectParams::default().plateau_fraction,
        }
    }

    /// The 5' mask length for a role.
    pub(crate) fn role_5p(&self, role: ReadRole) -> usize {
        match role {
            ReadRole::R1 => self.k_r1_5p,
            ReadRole::R2 => self.k_r2_5p,
            ReadRole::Se => self.k_se_5p,
        }
    }

    /// The detection threshold (`plateau_fraction × CpG plateau`) for a role's 5'
    /// end, or `None` when no plateau was estimable. Plots draw this as the cut
    /// line so the `0.9 × plateau` level is visible against the M-bias curve.
    pub(crate) fn threshold_5p(&self, role: ReadRole) -> Option<f64> {
        self.plateau_5p[role.index()].map(|pl| pl * self.plateau_fraction)
    }

    /// The `--mbias-plateau-fraction` the plan was learned with (for plot labels).
    pub(crate) fn plateau_fraction(&self) -> f64 {
        self.plateau_fraction
    }

    /// The single-end 3' mask length.
    pub(crate) fn k_se_3p(&self) -> usize {
        self.k_se_3p
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

/// Apply `plan` to every maskable record of one template, in place.
///
/// Maskable = has SEQ + real qualities and is not a secondary alignment, so
/// primary, supplementary, and unmapped records are all masked; secondaries are
/// skipped (see module docs). Three phases so cross-mate propagation can read
/// every record's own windows before any mutation: (1) each record's *own* mask
/// windows (5', plus SE 3' or the mate-5' mirror); (2) propagate masked reference
/// positions between the two mate sides by coverage, matched on contig; (3) write
/// the mask quality over all collected windows.
#[cfg(test)]
pub(crate) fn mask_template(plan: &MaskPlan, recs: &mut [RawRecord]) {
    let windows = compute_mask_windows(plan, recs);
    apply_windows(recs, &windows, plan.mask_quality);
}

/// Compute the stored-position mask windows for every record of one template,
/// without mutating anything (phases 1–2 of masking: each record's own windows,
/// then cross-mate reference-coverage propagation). `windows[i]` is the
/// set of half-open stored-position intervals to mask in `recs[i]`, empty for a
/// record that masks nothing (including secondaries).
///
/// Separated from the Q2 write ([`apply_windows`]) so the **same** geometry drives
/// both the masked output *and* the decision tally's exclusion — the reported
/// metrics then describe exactly the bases a downstream caller sees. The windows
/// can be interior (a propagated overhang), not just end trims, so callers must
/// treat them as a general interval set, never collapse them to a 5'/3' length.
pub(crate) fn compute_mask_windows(plan: &MaskPlan, recs: &[RawRecord]) -> Vec<MaskWindows> {
    let n = recs.len();

    // Does a mapped (non-secondary) record exist for each segment? A paired read
    // falls back to the mate-5' mirror when its mate has no mapped record to
    // propagate from.
    let (mut mapped_r1, mut mapped_r2) = (false, false);
    for rec in recs.iter() {
        let f = rec.flags();
        if has(f, FLAG_UNMAPPED | FLAG_SECONDARY) {
            continue;
        }
        mapped_r1 |= has(f, FLAG_FIRST_SEGMENT);
        mapped_r2 |= has(f, FLAG_LAST_SEGMENT);
    }

    // Phase 1 — each maskable record's own mask windows (stored positions).
    let mut windows: Vec<MaskWindows> = vec![MaskWindows::new(); n];
    for i in 0..n {
        if !maskable(&recs[i]) {
            continue;
        }
        let f = recs[i].flags();
        let seq_len = recs[i].l_seq() as usize;
        let reverse = has(f, FLAG_REVERSE);
        let role = read_role(f);
        let (left_hard, right_hard) = hard_clips(&recs[i]);

        // Own 5' (all roles), hard-clip-aware.
        if let Some(w) =
            five_prime_window(reverse, plan.role_5p(role), seq_len, left_hard, right_hard)
        {
            windows[i].push(w);
        }

        if role == ReadRole::Se {
            // Single-end: also the learned 3' end.
            if let Some(w) =
                three_prime_window(reverse, plan.k_se_3p, seq_len, left_hard, right_hard)
            {
                windows[i].push(w);
            }
        } else {
            // Paired read that can't rely on reference-coverage propagation
            // (phase 2) — its mate is unmapped/absent, or the read itself is
            // unmapped — mirrors the mate role's 5' length onto its 3' end, in
            // case it read through to the mate's 5'.
            let mate_present = if has(f, FLAG_FIRST_SEGMENT) { mapped_r2 } else { mapped_r1 };
            if has(f, FLAG_UNMAPPED) || has(f, FLAG_MATE_UNMAPPED) || !mate_present {
                let km = plan.role_5p(role.mate());
                if let Some(w) = three_prime_window(reverse, km, seq_len, left_hard, right_hard) {
                    windows[i].push(w);
                }
            }
        }
    }

    // Phase 2 — propagate masked reference positions between the two mate sides.
    // Collect each side's own-window reference ranges (carrying contig) before
    // applying anything, so propagated windows never feed back into the ranges.
    // A position masked on one side is then masked in any propagatable record of
    // the other side that covers it *on the same contig* — coverage- and
    // contig-driven, so it spans split-read supplementaries that cross an SV
    // breakpoint onto the mate's contig, not just same-contig primary mates.
    let mut r1_ranges: Vec<(i32, usize, usize)> = Vec::new();
    let mut r2_ranges: Vec<(i32, usize, usize)> = Vec::new();
    for i in 0..n {
        if windows[i].is_empty() || !propagatable(&recs[i]) {
            continue;
        }
        let tid = recs[i].ref_id();
        let side =
            if has(recs[i].flags(), FLAG_FIRST_SEGMENT) { &mut r1_ranges } else { &mut r2_ranges };
        for &(lo, hi) in &windows[i] {
            if let Some((ref_lo, ref_hi)) = ref_span_for_query_window(&recs[i], lo, hi) {
                side.push((tid, ref_lo, ref_hi));
            }
        }
    }
    for i in 0..n {
        if !propagatable(&recs[i]) {
            continue;
        }
        let tid = recs[i].ref_id();
        let other = if has(recs[i].flags(), FLAG_FIRST_SEGMENT) { &r2_ranges } else { &r1_ranges };
        for &(rid, ref_lo, ref_hi) in other {
            if rid == tid
                && let Some(w) = stored_window_for_ref(&recs[i], ref_lo, ref_hi)
            {
                windows[i].push(w);
            }
        }
    }

    windows
}

/// Write `mask_quality` over the precomputed mask windows (phase 3 of masking).
/// `windows[i]` applies to `recs[i]`; an empty entry leaves that record
/// untouched. Each interval is clamped to the stored quality length.
pub(crate) fn apply_windows(recs: &mut [RawRecord], windows: &[MaskWindows], mask_quality: u8) {
    for (i, w) in windows.iter().enumerate() {
        if w.is_empty() {
            continue;
        }
        let q = recs[i].quality_scores_mut();
        for &(lo, hi) in w {
            let hi = hi.min(q.len());
            for b in &mut q[lo.min(hi)..hi] {
                *b = mask_quality;
            }
        }
    }
}

/// Whether a record carries sequence and real base qualities — QUAL present, not
/// the `*` sentinel (stored as `0xFF`). The precondition for masking anything.
fn has_seq_and_qual(rec: &RawRecord) -> bool {
    rec.l_seq() as usize > 0 && rec.quality_scores().first().is_some_and(|&q| q != 0xFF)
}

/// Whether a record is eligible for masking: it has SEQ + real qualities and is
/// not a secondary alignment.
fn maskable(rec: &RawRecord) -> bool {
    !has(rec.flags(), FLAG_SECONDARY) && has_seq_and_qual(rec)
}

/// Whether a record can take part in cross-mate reference-coverage propagation:
/// maskable, mapped (so it has reference coordinates), and paired.
fn propagatable(rec: &RawRecord) -> bool {
    let f = rec.flags();
    maskable(rec) && !has(f, FLAG_UNMAPPED) && has(f, FLAG_PAIRED) && rec.ref_id() >= 0
}

/// Leading and trailing hard-clip lengths (CIGAR `H` ops). Hard-clipped bases are
/// absent from SEQ, so they offset the read's sequencing-cycle frame from the
/// stored-position frame. `H` is only ever the first and/or last op; a single-op
/// CIGAR can't be both, so trailing is taken only when there are at least two ops.
fn hard_clips(rec: &RawRecord) -> (usize, usize) {
    let mut ops = rec.cigar_ops_iter();
    let Some(first) = ops.next() else {
        return (0, 0);
    };
    let leading = if first & 0xf == 5 { (first >> 4) as usize } else { 0 };
    let mut last = first;
    let mut count = 1usize;
    for op in ops {
        last = op;
        count += 1;
    }
    let trailing = if count > 1 && last & 0xf == 5 { (last >> 4) as usize } else { 0 };
    (leading, trailing)
}

/// Stored half-open window covering the `k` sequencing cycles nearest the read's
/// **5' end**, accounting for any 5' cycles hard-clipped out of SEQ. `None` when
/// `k == 0` or every targeted cycle was hard-clipped away.
///
/// SEQ is stored genomic-forward, so a forward read's 5' is the low stored end
/// (clipped by `left_hard`) and a reverse read's is the high stored end (clipped
/// by `right_hard`). With no hard clips this is just `[0, k)` / `[len - k, len)`.
fn five_prime_window(
    reverse: bool,
    k: usize,
    len: usize,
    left_hard: usize,
    right_hard: usize,
) -> Option<(usize, usize)> {
    if k == 0 {
        return None;
    }
    if reverse {
        let start = (len + right_hard).saturating_sub(k);
        (start < len).then_some((start, len))
    } else {
        let end = k.saturating_sub(left_hard).min(len);
        (end > 0).then_some((0, end))
    }
}

/// Stored window for the `k` cycles nearest the read's **3' end** — the mirror of
/// [`five_prime_window`], since the 3' end is the 5' end of the opposite
/// orientation.
fn three_prime_window(
    reverse: bool,
    k: usize,
    len: usize,
    left_hard: usize,
    right_hard: usize,
) -> Option<(usize, usize)> {
    five_prime_window(!reverse, k, len, left_hard, right_hard)
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
    use crate::record::{FLAG_LAST_SEGMENT, FLAG_SUPPLEMENTARY};
    use crate::sam_reader::SamReader;

    /// Read a complete SAM document (header + records) into `RawRecord`s.
    fn read_sam(sam: String) -> Vec<RawRecord> {
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

    /// Parse SAM lines against a single `chr1` contig.
    fn parse(records: &[&str], contig_len: usize) -> Vec<RawRecord> {
        parse_multi(&[("chr1", contig_len)], records)
    }

    /// Parse SAM lines against an explicit set of `(name, length)` contigs.
    fn parse_multi(contigs: &[(&str, usize)], records: &[&str]) -> Vec<RawRecord> {
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for (name, len) in contigs {
            sam.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
        }
        for r in records {
            sam.push_str(r);
            sam.push('\n');
        }
        read_sam(sam)
    }

    fn line(qname: &str, flag: u16, pos: u32, cigar: &str, seq: &str) -> String {
        line_on(qname, flag, "chr1", pos, cigar, seq)
    }

    /// A mapped SAM line on an explicit contig (`I`-quality, MAPQ 60).
    fn line_on(qname: &str, flag: u16, rname: &str, pos: u32, cigar: &str, seq: &str) -> String {
        let qual = "I".repeat(seq.len());
        format!("{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
    }

    /// An unmapped SAM line (RNAME `*`, POS 0, CIGAR `*`, MAPQ 0).
    fn line_unmapped(qname: &str, flag: u16, seq: &str) -> String {
        let qual = "I".repeat(seq.len());
        format!("{qname}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
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
    /// propagation occurs — only each read's own 5' is masked.
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

    /// A supplementary alignment is masked too, but its 5' hard-clipped cycles are
    /// absent from SEQ, so the stored window shifts by the leading hard-clip
    /// length: a K=8 5' mask with 5 cycles hard-clipped masks only the 3 present.
    #[test]
    fn hard_clipped_supplementary_shifts_5p_window() {
        let mut recs = parse(&[&line("s", R1 | FLAG_SUPPLEMENTARY, 1, "5H10M", "CACACACACA")], 40);
        mask_template(&MaskPlan::explicit(8, 0, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[0..3], &[2, 2, 2], "8-cycle 5' mask − 5 hard-clipped = 3 stored");
        assert!(q[3..].iter().all(|&b| b == 40), "rest untouched");
    }

    /// When the whole 5' mask falls inside the hard-clipped region, nothing in
    /// this record's SEQ is masked.
    #[test]
    fn hard_clipped_supplementary_fully_past_window_masks_nothing() {
        let mut recs = parse(&[&line("s", R1 | FLAG_SUPPLEMENTARY, 1, "10H10M", "CACACACACA")], 40);
        mask_template(&MaskPlan::explicit(8, 0, 2), &mut recs);
        assert!(quals(&recs[0]).iter().all(|&b| b == 40), "5' mask entirely hard-clipped away");
    }

    /// An unmapped paired read is masked defensively at both ends: own 5' by its
    /// role's length, and the 3' by the mate role's 5' length (mirror), since it
    /// can't recover the mate's masked positions from coverage.
    #[test]
    fn unmapped_paired_read_masks_own_5p_and_mate_mirror_3p() {
        let mut recs = parse(&[&line_unmapped("u", R1 | FLAG_UNMAPPED, "CACACACACA")], 20);
        mask_template(&MaskPlan::explicit(3, 4, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[0..3], &[2, 2, 2], "own 5' (k_r1)");
        assert!(q[3..6].iter().all(|&b| b == 40), "middle kept");
        assert_eq!(&q[6..10], &[2, 2, 2, 2], "mate-5' (k_r2) mirror on 3'");
    }

    /// Reverse unmapped read: 5' is the high stored end, 3' the low end (we trust
    /// the reverse flag even with no alignment).
    #[test]
    fn unmapped_reverse_read_masks_correct_ends() {
        let mut recs =
            parse(&[&line_unmapped("u", R1 | FLAG_UNMAPPED | FLAG_REVERSE, "CACACACACA")], 20);
        mask_template(&MaskPlan::explicit(3, 4, 2), &mut recs);
        let q = quals(&recs[0]);
        assert_eq!(&q[7..10], &[2, 2, 2], "own 5' at high stored end");
        assert_eq!(&q[0..4], &[2, 2, 2, 2], "mate-5' mirror on low (3') end");
        assert!(q[4..7].iter().all(|&b| b == 40), "middle kept");
    }

    /// Secondary alignments are never masked, even with SEQ and qualities present.
    #[test]
    fn secondary_alignment_is_never_masked() {
        let mut recs = parse(&[&line("x", R1 | FLAG_SECONDARY, 1, "10M", "CACACACACA")], 20);
        mask_template(&MaskPlan::explicit(5, 5, 2), &mut recs);
        assert!(quals(&recs[0]).iter().all(|&b| b == 40), "secondary untouched");
    }

    /// A record with no base qualities (QUAL = `*`) is skipped — we never
    /// fabricate qualities, so the missing-quality sentinel is left intact.
    #[test]
    fn record_without_qualities_is_skipped() {
        let mut recs = parse(&["n\t65\tchr1\t1\t60\t10M\t*\t0\t0\tCACACACACA\t*"], 20);
        mask_template(&MaskPlan::explicit(5, 0, 2), &mut recs);
        assert!(quals(&recs[0]).iter().all(|&b| b == 0xFF), "missing QUAL left as sentinel");
    }

    /// Split read across an SV breakpoint: R1's primary maps to chr1 and R2's to
    /// chr3, but each read's *supplementary* lands on the mate's contig. Reference
    /// coverage propagation must cross contigs so each 5'-hard-clipped
    /// supplementary still gets masked at the mate's biased 5' positions.
    #[test]
    fn sv_split_read_propagates_across_contigs() {
        let r1_prim =
            line_on("p", FLAG_PAIRED | FLAG_FIRST_SEGMENT, "chr1", 1, "10M", "CACACACACA");
        let r1_supp = line_on(
            "p",
            FLAG_PAIRED | FLAG_FIRST_SEGMENT | FLAG_SUPPLEMENTARY,
            "chr3",
            1,
            "10H10M",
            "CACACACACA",
        );
        let r2_prim = line_on("p", FLAG_PAIRED | FLAG_LAST_SEGMENT, "chr3", 1, "10M", "CACACACACA");
        let r2_supp = line_on(
            "p",
            FLAG_PAIRED | FLAG_LAST_SEGMENT | FLAG_SUPPLEMENTARY,
            "chr1",
            1,
            "10H10M",
            "CACACACACA",
        );
        let mut recs =
            parse_multi(&[("chr1", 40), ("chr3", 40)], &[&r1_prim, &r1_supp, &r2_prim, &r2_supp]);
        mask_template(&MaskPlan::explicit(3, 4, 2), &mut recs); // k_r1=3, k_r2=4

        // R1 primary (chr1): own 5' only; R2's ranges are on chr3, so no propagation.
        assert_eq!(&quals(&recs[0])[0..3], &[2, 2, 2], "R1 primary own 5'");
        assert!(quals(&recs[0])[3..].iter().all(|&b| b == 40));
        // R1 supplementary (chr3): 5'-hard-clipped (no own window) → masked only by
        // R2's chr3 5' positions [0,4) via cross-contig propagation.
        assert_eq!(&quals(&recs[1])[0..4], &[2, 2, 2, 2], "R1 supp masked by R2's chr3 5'");
        assert!(quals(&recs[1])[4..].iter().all(|&b| b == 40));
        // R2 primary (chr3): own 5' only.
        assert_eq!(&quals(&recs[2])[0..4], &[2, 2, 2, 2], "R2 primary own 5'");
        assert!(quals(&recs[2])[4..].iter().all(|&b| b == 40));
        // R2 supplementary (chr1): masked by R1's chr1 5' positions [0,3).
        assert_eq!(&quals(&recs[3])[0..3], &[2, 2, 2], "R2 supp masked by R1's chr1 5'");
        assert!(quals(&recs[3])[3..].iter().all(|&b| b == 40));
    }
}
