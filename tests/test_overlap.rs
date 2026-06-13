//! PE-overlap deduplication: reference positions covered by both mates of an
//! overlapping proper pair are counted once — the overlap is split at its
//! midpoint and each mate keeps the half nearer its own 5' end (higher base
//! quality). This matters most for short-insert / cfDNA libraries, where the
//! overlap can be most of the read.

mod helpers;
use helpers::*;

// CpA reference: top-strand C at even offsets.
const REF20: &str = "CACACACACACACACACACA"; // 20 bp, C at 0,2,..,18
const REF10: &str = "CACACACACA"; // 10 bp, C at 0,2,4,6,8

fn ot_r1() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT | FLAG_MATE_REVERSE
}
fn ot_r2() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_LAST_SEGMENT | FLAG_REVERSE
}
fn ca_total(env: &TestEnv) -> u64 {
    genome_stats(&env.stats)["CA_total"].parse().unwrap()
}
fn ca_unconv(env: &TestEnv) -> u64 {
    genome_stats(&env.stats)["CA_unconv"].parse().unwrap()
}

#[test]
fn overlap_dedup_keeps_template_below_threshold() {
    // Fully-overlapping mates (both at pos 1, 6M over ref[0,6)), each with 2
    // unconverted C's at the SAME positions (C@0,2; C@4 converted). Counting the
    // overlap once leaves 2 unique unconverted C's (< 3) → NOT tagged. Without
    // dedup the same data would naively double to 4 and cross the threshold, so
    // this confirms dedup is what holds the decision.
    let reference = RefBuilder::new().contig("chr1", REF10);
    let sam = SamBuilder::new()
        .sq("chr1", REF10.len())
        .record("p", ot_r1(), "chr1", 1, "6M", "CACATA", &q40(6))
        .record("p", ot_r2(), "chr1", 1, "6M", "CACATA", &q40(6));

    let env = TestEnv::new();
    let deduped = run_ok(&sam, &reference, &env, &[]);
    assert!(
        deduped.iter().all(|r| !has_tag(r, [b'X', b'X'])),
        "with overlap dedup, 2 unique unconverted < 3 → not tagged"
    );
}

#[test]
fn each_overlapped_reference_base_counted_once() {
    // R1 covers ref[0,10) (C@0,2,4,6,8); R2 covers ref[5,15) (C@6,8,10,12,14).
    // Overlap ref[5,10) holds C@6,8 and is split at its midpoint between the
    // mates. Deduped total CpA sites = 8 (naive double-counting would give 10).
    let reference = RefBuilder::new().contig("chr1", REF20);
    let sam = SamBuilder::new()
        .sq("chr1", REF20.len())
        .record("p", ot_r1(), "chr1", 1, "10M", "CACACACACA", &q40(10)) // ref[0,10)
        .record("p", ot_r2(), "chr1", 6, "10M", "ACACACACAC", &q40(10)); // ref[5,15)

    let env = TestEnv::new();
    let stats = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert_eq!(ca_total(&env), 8, "overlapped C@6,8 counted once");
}

#[test]
fn overlap_split_assigns_each_half_to_the_nearer_mate() {
    // Fully-overlapping mates over ref[0,10). R1 (forward) reads every C as
    // unconverted; R2 (reverse) reads every C as converted. The overlap splits at
    // the midpoint (5): R1 keeps [0,5) → C@0,2,4 (3 unconverted), R2 keeps [5,10)
    // → C@6,8 (converted). So 5 distinct CpA sites, 3 unconverted — whereas the
    // old "R1 takes the whole overlap" would have scored all 5 unconverted.
    let reference = RefBuilder::new().contig("chr1", REF10);
    let sam = SamBuilder::new()
        .sq("chr1", REF10.len())
        .record("p", ot_r1(), "chr1", 1, "10M", "CACACACACA", &q40(10))
        .record("p", ot_r2(), "chr1", 1, "10M", "TATATATATA", &q40(10));

    let env = TestEnv::new();
    let stats = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert_eq!(ca_total(&env), 5, "5 distinct CpA sites across the split overlap");
    assert_eq!(ca_unconv(&env), 3, "only R1's half [0,5) contributes unconverted C's");
}

#[test]
fn non_overlapping_pair_is_unchanged() {
    // R1 ref[0,10), R2 ref[10,20): adjacent, no overlap. Both modes tally all 10.
    let reference = RefBuilder::new().contig("chr1", REF20);
    let sam = SamBuilder::new()
        .sq("chr1", REF20.len())
        .record("p", ot_r1(), "chr1", 1, "10M", "CACACACACA", &q40(10)) // ref[0,10)
        .record("p", ot_r2(), "chr1", 11, "10M", "CACACACACA", &q40(10)); // ref[10,20)

    let env = TestEnv::new();
    let stats = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert_eq!(ca_total(&env), 10, "non-overlapping mates both fully counted");
}

#[test]
fn single_end_reads_are_unaffected_by_dedup() {
    // Two independent SE reads at the same locus are NOT mates — both counted.
    let reference = RefBuilder::new().contig("chr1", REF10);
    let sam = SamBuilder::new()
        .sq("chr1", REF10.len())
        .record("a", 0, "chr1", 1, "10M", "CACACACACA", &q40(10))
        .record("b", 0, "chr1", 1, "10M", "CACACACACA", &q40(10));

    let env = TestEnv::new();
    let stats = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    // Two separate templates, 5 CpA each → 10 total; dedup never applies to SE.
    assert_eq!(ca_total(&env), 10);
}
