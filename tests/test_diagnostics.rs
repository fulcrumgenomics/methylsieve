//! Diagnostic counters and pass-through behavior: unmapped templates, zero-site
//! templates, and the count-threshold boundary.

mod helpers;
use helpers::*;

const REF: &str = "CACACACACA"; // 5 CpA C's at 0,2,4,6,8

#[test]
fn unmapped_template_passes_through_untouched() {
    // A read pair with both ends unmapped: no tag, no QC-fail, counted as
    // unmapped, never evaluated.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let stats = env.stats.to_str().unwrap().to_string();
    let r1 = FLAG_PAIRED | FLAG_UNMAPPED | FLAG_MATE_UNMAPPED | FLAG_FIRST_SEGMENT;
    let r2 = FLAG_PAIRED | FLAG_UNMAPPED | FLAG_MATE_UNMAPPED | FLAG_LAST_SEGMENT;
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("u", r1, "*", 0, "*", "CACACACACA", &q40(10))
        .record("u", r2, "*", 0, "*", "CACACACACA", &q40(10));

    let recs = run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert_eq!(recs.len(), 2);
    for rec in &recs {
        assert!(!has_tag(rec, [b'X', b'X']), "unmapped reads must not be tagged");
        assert_eq!(u16::from(rec.flags()) & FLAG_QC_FAIL, 0);
    }
    let g = genome_stats(&env.stats);
    assert_eq!(g["unmapped_templates"], "1");
    assert_eq!(g["n_evaluated"], "0");
}

#[test]
fn zero_site_template_is_counted_and_not_tagged() {
    // A read over an all-A region has no monitored C — zero sites.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", "AAAAAAAAAA");
    let stats = env.stats.to_str().unwrap().to_string();
    let sam =
        SamBuilder::new().sq("chr1", 10).record("z", 0, "chr1", 1, "10M", "AAAAAAAAAA", &q40(10));

    let recs = run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert!(!has_tag(&recs[0], [b'X', b'X']));
    let g = genome_stats(&env.stats);
    assert_eq!(g["zero_site_templates"], "1");
    assert_eq!(g["n_evaluated"], "0");
}

#[test]
fn count_threshold_boundary_is_inclusive() {
    // Exactly 2 unconverted C's: tagged with --max-unconverted-count 2, not with 3.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACATATATA",
        &q40(10), // C@0,2 unconverted
    );

    let tagged_at_2 = run_ok(&sam, &reference, &env, &["--max-unconverted-count", "2"]);
    assert!(has_tag(&tagged_at_2[0], [b'X', b'X']), "2 ≥ 2 → tagged");

    let env2 = TestEnv::new();
    let not_at_3 = run_ok(&sam, &reference, &env2, &["--max-unconverted-count", "3"]);
    assert!(!has_tag(&not_at_3[0], [b'X', b'X']), "2 < 3 → not tagged");
}

#[test]
fn count_of_zero_is_rejected() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam =
        SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "10M", REF, &q40(10));
    let r = run_methylsieve(&sam, &reference, &env, &["--max-unconverted-count", "0"]);
    assert!(!r.status_ok, "count 0 should be rejected");
    assert!(r.stderr.contains("max-unconverted-count"), "got: {}", r.stderr);
}

#[test]
fn deletion_skips_reference_cytosines() {
    // CIGAR 4M2D4M over ref "CACACACACA": the 2D skips ref positions 4,5 (a C at
    // 4). The read's 8 bases align to ref 0-3 and 6-9. Monitored C's: ref 0,2
    // (first block) and 6,8 (second block) = 4; the deleted C@4 is NOT tallied.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let stats = env.stats.to_str().unwrap().to_string();
    // Read bases (8) for ref 0-3 ("CACA") then ref 6-9 ("CACA"), all matching.
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "d",
        0,
        "chr1",
        1,
        "4M2D4M",
        "CACACACA",
        &q40(8),
    );
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    let g = genome_stats(&env.stats);
    assert_eq!(g["CA_total"], "4", "deleted C@4 must not be tallied");
    assert_eq!(g["CA_unconv"], "4");
}
