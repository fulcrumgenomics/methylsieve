//! Per-record strand determination (`monitor_C = (R1 or unpaired) XOR reverse`).
//!
//! These tests are designed to FAIL under the naive alternatives:
//! * "use each record's own 0x10" (would mis-bucket R2), and
//! * "propagate primary R1's orientation to every record" (would mis-bucket a
//!   strand-flipping supplementary).

mod helpers;
use helpers::*;

// Top-strand region: C at 0,2,4,6,8 — each CpA. No G anywhere.
const TOP_REF: &str = "CACACACACA";
// Bottom-strand region: G at 1,3,5,7,9 — each preceded by T, so CpA-equivalent.
// No C anywhere.
const BOT_REF: &str = "TGTGTGTGTG";

// 20 bp variants so the two mates can sit at non-overlapping loci (R1 ref[0,10),
// R2 ref[10,20)) — otherwise PE-overlap dedup would collapse the duplicated
// evidence and defeat the strand discriminator.
const TOP_REF20: &str = "CACACACACACACACACACA";
const BOT_REF20: &str = "TGTGTGTGTGTGTGTGTGTG";

#[test]
fn ot_pair_both_reads_monitor_ref_c() {
    // OT pair: R1 forward, R2 reverse. Both must monitor ref C. R1 (ref[0,10))
    // contributes 2 unconverted C's and R2 (ref[10,20)) another 2 → 4 ≥ 3 →
    // tagged. Under "R2 uses its own reverse bit" R2 would monitor G (none
    // present) → only 2 → not tagged.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", TOP_REF20);
    let r1_flag = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT | FLAG_MATE_REVERSE;
    let r2_flag = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_LAST_SEGMENT | FLAG_REVERSE;
    // SEQ (ref-forward orientation): C at the first two even positions, T else.
    let two_unconv = "CACATATATA";
    let sam = SamBuilder::new()
        .sq("chr1", TOP_REF20.len())
        .record("p", r1_flag, "chr1", 1, "10M", two_unconv, &q40(10)) // ref[0,10)
        .record("p", r2_flag, "chr1", 11, "10M", two_unconv, &q40(10)); // ref[10,20)

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 2);
    for rec in &recs {
        assert!(has_tag(rec, *b"XX"), "both records of OT pair should be tagged");
    }
}

#[test]
fn ob_pair_both_reads_monitor_ref_g() {
    // OB pair: R1 reverse, R2 forward. Both must monitor ref G. Mates at
    // non-overlapping loci so their evidence sums (2 + 2 ≥ 3).
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", BOT_REF20);
    let r1_flag = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT | FLAG_REVERSE;
    let r2_flag = FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_LAST_SEGMENT | FLAG_MATE_REVERSE;
    // SEQ (ref-forward orientation): G at the first two odd positions, A else.
    let two_unconv = "TGTGTATATA";
    let sam = SamBuilder::new()
        .sq("chr1", BOT_REF20.len())
        .record("p", r1_flag, "chr1", 1, "10M", two_unconv, &q40(10)) // ref[0,10)
        .record("p", r2_flag, "chr1", 11, "10M", two_unconv, &q40(10)); // ref[10,20)

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 2);
    for rec in &recs {
        assert!(has_tag(rec, *b"XX"), "both records of OB pair should be tagged");
    }
}

#[test]
fn reverse_supplementary_flips_to_ref_g() {
    // Single-end molecule: a forward primary over a C region (1 unconverted C)
    // plus a reverse-mapped supplementary over a G region (2 unconverted G).
    // Per-record strand: primary→monitor C, supplementary→monitor G, summing to
    // 3 ≥ 3 → tagged. Under "all records inherit primary R1's forward
    // orientation" the supplementary would monitor C (no C in the G region) →
    // only 1 → NOT tagged. So this asserts the per-record fix.
    let env = TestEnv::new();
    // chr1[0..10] = C region (CpA C's at 0,2,4,6,8); chr1[10..20] = G region.
    let combined = format!("{TOP_REF}{BOT_REF}");
    let reference = RefBuilder::new().contig("chr1", &combined);

    let primary_flag = 0u16; // unpaired, forward
    let supp_flag = FLAG_SUPPLEMENTARY | FLAG_REVERSE; // unpaired, reverse supplementary
    let sam = SamBuilder::new()
        .sq("chr1", combined.len())
        // primary: 1 unconverted C (only position 0 reads C).
        .record("t1", primary_flag, "chr1", 1, "10M", "CATATATATA", &q40(10))
        // supplementary at pos 11 (G region): 2 unconverted G (positions 1,3).
        .record("t1", supp_flag, "chr1", 11, "10M", "TGTGTATATA", &q40(10));

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 2);
    for rec in &recs {
        assert!(
            has_tag(rec, *b"XX"),
            "tag must propagate to all records (incl. the supplementary)"
        );
    }
}

#[test]
fn ignore_supplementary_evidence_excludes_it_from_decision() {
    // Same setup as above, but suppress supplementary evidence: now only the
    // primary's 1 unconverted C counts → below threshold → not tagged.
    let env = TestEnv::new();
    let combined = format!("{TOP_REF}{BOT_REF}");
    let reference = RefBuilder::new().contig("chr1", &combined);
    let sam = SamBuilder::new()
        .sq("chr1", combined.len())
        .record("t1", 0, "chr1", 1, "10M", "CATATATATA", &q40(10))
        .record("t1", FLAG_SUPPLEMENTARY | FLAG_REVERSE, "chr1", 11, "10M", "TGTGTATATA", &q40(10));

    let recs = run_ok(&sam, &reference, &env, &["--ignore-supplementary-evidence"]);
    for rec in &recs {
        assert!(!has_tag(rec, *b"XX"), "no record should be tagged when below threshold");
    }
}
