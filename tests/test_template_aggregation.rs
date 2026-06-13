//! Evidence is aggregated across all evaluated records of a template, and the
//! decision is propagated to every record; secondary alignments never
//! contribute evidence.

mod helpers;
use helpers::*;

// 20 bp so R1 (ref[0,10)) and R2 (ref[10,20)) sit at non-overlapping loci and
// their evidence sums without PE-overlap dedup collapsing it.
const REF: &str = "CACACACACACACACACACA";

fn ot_r1() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT | FLAG_MATE_REVERSE
}
fn ot_r2() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_LAST_SEGMENT | FLAG_REVERSE
}

#[test]
fn aggregates_r1_r2_and_supplementary_then_tags_all() {
    // 1 (R1) + 2 (R2) + 0 (supp) = 3 ≥ 3 → unconverted; tag must land on all
    // three records, supplementary included.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("t", ot_r1(), "chr1", 1, "10M", "CATATATATA", &q40(10)) // ref[0,10): 1 unconv
        .record("t", ot_r2(), "chr1", 11, "10M", "CACATATATA", &q40(10)) // ref[10,20): 2 unconv
        .record(
            "t",
            FLAG_PAIRED | FLAG_FIRST_SEGMENT | FLAG_SUPPLEMENTARY,
            "chr1",
            1,
            "10M",
            "TATATATATA",
            &q40(10),
        ); // 0 unconv

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 3);
    for rec in &recs {
        assert!(has_tag(rec, [b'X', b'X']), "tag must propagate to every record");
    }
}

#[test]
fn secondary_alignments_do_not_contribute_evidence() {
    // R1 (1) + R2 (1) = 2 < 3. A secondary record carrying 5 unconverted C's
    // must NOT push the template over the threshold.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("t", ot_r1(), "chr1", 1, "10M", "CATATATATA", &q40(10)) // ref[0,10): 1 unconv
        .record("t", ot_r2(), "chr1", 11, "10M", "CATATATATA", &q40(10)) // ref[10,20): 1 unconv
        .record(
            "t",
            FLAG_PAIRED | FLAG_FIRST_SEGMENT | FLAG_SECONDARY,
            "chr1",
            1,
            "10M",
            "CACACACACA",
            &q40(10),
        ); // 5 unconv but secondary → ignored

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 3);
    for rec in &recs {
        assert!(!has_tag(rec, [b'X', b'X']), "secondary evidence must not trigger a tag");
    }
}
