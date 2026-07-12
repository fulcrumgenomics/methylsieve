//! Low-level BAM record geometry shared by the tally engine ([`crate::sieve`])
//! and the masker ([`crate::mask`]): SAM flag constants, flag predicates, the
//! per-record monitored-strand / read-role rules, and CIGAR-aware mapping
//! between stored read positions and reference coordinates.
//!
//! These are per-*record* helpers (called O(records), not O(bases)), so they sit
//! outside the per-base tally hot loop in `sieve`.

use fgumi_raw_bam::RawRecord;

use crate::mbias::ReadRole;

// ── SAM flag constants ──────────────────────────────────────────────────────

/// SAM FLAG 0x1: template has multiple segments (paired-end).
pub(crate) const FLAG_PAIRED: u16 = 0x1;
/// SAM FLAG 0x4: segment unmapped.
pub(crate) const FLAG_UNMAPPED: u16 = 0x4;
/// SAM FLAG 0x8: the mate / next segment is unmapped.
pub(crate) const FLAG_MATE_UNMAPPED: u16 = 0x8;
/// SAM FLAG 0x10: this segment is on the reverse strand.
pub(crate) const FLAG_REVERSE: u16 = 0x10;
/// SAM FLAG 0x40: first segment in the template (R1).
pub(crate) const FLAG_FIRST_SEGMENT: u16 = 0x40;
/// SAM FLAG 0x80: last segment in the template (R2).
pub(crate) const FLAG_LAST_SEGMENT: u16 = 0x80;
/// SAM FLAG 0x100: secondary alignment.
pub(crate) const FLAG_SECONDARY: u16 = 0x100;
/// SAM FLAG 0x200: QC-fail (the bit we OR in for unconverted templates).
pub(crate) const FLAG_QC_FAIL: u16 = 0x200;
/// SAM FLAG 0x800: supplementary (chimeric) alignment.
pub(crate) const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// Whether `flags` has **any** bit in `bit` set (`flags & bit != 0`). For a
/// single-bit `bit` this is the obvious "is this flag set?"; for an OR'd mask it
/// is "any of these flags set" — e.g. `has(f, FLAG_SECONDARY | FLAG_SUPPLEMENTARY)`
/// tests secondary-or-supplementary.
#[inline]
pub(crate) fn has(flags: u16, bit: u16) -> bool {
    flags & bit != 0
}

/// Whether a record monitors reference C (top, `true`) or G (bottom, `false`),
/// per the per-record MethylDackel rule: treat single-end and R1 the same, then
/// XOR with the record's own reverse bit.
#[inline]
pub(crate) fn monitor_c_of(flags: u16) -> bool {
    let treat_as_read1 = has(flags, FLAG_FIRST_SEGMENT) || !has(flags, FLAG_PAIRED);
    treat_as_read1 ^ has(flags, FLAG_REVERSE)
}

/// Read role for M-bias bucketing: unpaired reads are [`ReadRole::Se`] (their
/// far template end is unknown, so both ends are analyzed); paired reads are
/// R1/R2 by segment bit. Distinct from [`monitor_c_of`], which folds single-end
/// into R1.
#[inline]
pub(crate) fn read_role(flags: u16) -> ReadRole {
    if !has(flags, FLAG_PAIRED) {
        ReadRole::Se
    } else if has(flags, FLAG_LAST_SEGMENT) {
        ReadRole::R2
    } else {
        ReadRole::R1
    }
}

/// Whether a record is a primary mapped alignment: mapped, not secondary, not
/// supplementary. This is the M-bias accumulation population (masking targets a
/// broader set — see `maskable` in [`crate::mask`]).
#[inline]
pub(crate) fn is_primary_mapped(flags: u16) -> bool {
    !has(flags, FLAG_UNMAPPED | FLAG_SECONDARY | FLAG_SUPPLEMENTARY)
}

/// Reference half-open interval `[start, end)` covered by the alignment of the
/// stored query positions in `[q_lo, q_hi)`. Returns `None` if no aligned base
/// falls in that window (e.g. a fully soft-clipped end). Insertions/soft-clips
/// inside the window consume query budget but contribute no reference positions;
/// a deletion inside the window stretches the returned span across the gap, so
/// the result is a single contiguous interval. Assumes `rec.pos() >= 0`.
pub(crate) fn ref_span_for_query_window(
    rec: &RawRecord,
    q_lo: usize,
    q_hi: usize,
) -> Option<(usize, usize)> {
    if q_lo >= q_hi {
        return None;
    }
    let mut qpos = 0usize;
    let mut rpos = rec.pos().max(0) as usize;
    let mut lo = usize::MAX;
    let mut hi = 0usize;
    for op in rec.cigar_ops_iter() {
        let len = (op >> 4) as usize;
        match op & 0xf {
            // M, =, X — aligned 1:1; intersect this op's query range with the window.
            0 | 7 | 8 => {
                let a = qpos.max(q_lo);
                let b = (qpos + len).min(q_hi);
                if a < b {
                    lo = lo.min(rpos + (a - qpos));
                    hi = hi.max(rpos + (b - qpos));
                }
                qpos += len;
                rpos += len;
            }
            // I, S — consume query only.
            1 | 4 => qpos += len,
            // D, N — consume reference only.
            2 | 3 => rpos += len,
            // H, P — consume neither.
            _ => {}
        }
        if qpos >= q_hi {
            break;
        }
    }
    (lo < hi).then_some((lo, hi))
}

/// Reference span of the `n` sequenced bases at the read's 5' end (sequencing
/// order). For a forward record the 5' end is the low stored-position end; for a
/// reverse record SEQ is stored reverse-complemented, so the 5' end is the high
/// stored-position end. Returns `None` if `n == 0` or those bases include no
/// aligned position.
pub(crate) fn five_prime_ref_span(rec: &RawRecord, n: usize) -> Option<(usize, usize)> {
    if n == 0 {
        return None;
    }
    let seq_len = rec.l_seq() as usize;
    let (q_lo, q_hi) =
        if has(rec.flags(), FLAG_REVERSE) { (seq_len.saturating_sub(n), seq_len) } else { (0, n) };
    ref_span_for_query_window(rec, q_lo, q_hi)
}

/// Reference span of the `n` sequenced bases at the read's 3' end — the opposite
/// end from [`five_prime_ref_span`]. Used only for a lone (single-end or orphan)
/// read, whose far template terminus can't be located, so both ends are trimmed.
pub(crate) fn three_prime_ref_span(rec: &RawRecord, n: usize) -> Option<(usize, usize)> {
    if n == 0 {
        return None;
    }
    let seq_len = rec.l_seq() as usize;
    let (q_lo, q_hi) =
        if has(rec.flags(), FLAG_REVERSE) { (0, n) } else { (seq_len.saturating_sub(n), seq_len) };
    ref_span_for_query_window(rec, q_lo, q_hi)
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, Cursor};

    use super::*;
    use crate::sam_reader::SamReader;

    /// Parse one-contig SAM lines into `RawRecord`s via methylsieve's SAM reader.
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

    fn line(flag: u16, pos: u32, cigar: &str, seq: &str) -> String {
        let q = "I".repeat(seq.len());
        format!("q\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{q}")
    }

    #[test]
    fn monitor_strand_rule() {
        // R1 fwd → C; R1 rev → G; R2 fwd → G; R2 rev → C; SE fwd → C.
        assert!(monitor_c_of(FLAG_PAIRED | FLAG_FIRST_SEGMENT));
        assert!(!monitor_c_of(FLAG_PAIRED | FLAG_FIRST_SEGMENT | FLAG_REVERSE));
        assert!(!monitor_c_of(FLAG_PAIRED | FLAG_LAST_SEGMENT));
        assert!(monitor_c_of(FLAG_PAIRED | FLAG_LAST_SEGMENT | FLAG_REVERSE));
        assert!(monitor_c_of(0));
    }

    #[test]
    fn read_role_rule() {
        assert_eq!(read_role(0), ReadRole::Se);
        assert_eq!(read_role(FLAG_PAIRED | FLAG_FIRST_SEGMENT), ReadRole::R1);
        assert_eq!(read_role(FLAG_PAIRED | FLAG_LAST_SEGMENT), ReadRole::R2);
    }

    #[test]
    fn five_prime_span_forward_and_reverse() {
        let fwd = parse(&[&line(0, 1, "10M", "CACACACACA")], 30);
        assert_eq!(five_prime_ref_span(&fwd[0], 3), Some((0, 3)));
        assert_eq!(three_prime_ref_span(&fwd[0], 3), Some((7, 10)));

        // Reverse: SEQ is stored forward-genomic, so the 5' end is the HIGH end.
        let rev = parse(&[&line(FLAG_REVERSE, 1, "10M", "CACACACACA")], 30);
        assert_eq!(five_prime_ref_span(&rev[0], 3), Some((7, 10)));
        assert_eq!(three_prime_ref_span(&rev[0], 3), Some((0, 3)));
    }

    #[test]
    fn five_prime_span_skips_leading_soft_clip() {
        let recs = parse(&[&line(0, 1, "3S7M", "GGGCACACAC")], 30);
        assert_eq!(five_prime_ref_span(&recs[0], 3), None);
        assert_eq!(five_prime_ref_span(&recs[0], 5), Some((0, 2)));
    }

    #[test]
    fn ref_span_stretches_across_deletion() {
        let recs = parse(&[&line(0, 1, "2M3D4M", "CACGTA")], 30);
        assert_eq!(five_prime_ref_span(&recs[0], 4), Some((0, 7)));
    }
}
