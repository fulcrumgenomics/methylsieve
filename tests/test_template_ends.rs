//! `--ignore-template-ends` for paired-end reads: trimming is by genomic
//! position at the fragment termini (the 5' sequenced ends of R1 and R2), not by
//! read position. Verified through the per-context tallies in the `--stats` TSV.

mod helpers;
use helpers::*;

// CpA reference: top-strand C at even offsets. In this synthetic FR layout both
// mates monitor reference C (MethylDackel getStrand), so a C-only reference
// exercises both reads.
fn ot_r1() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT | FLAG_MATE_REVERSE
}
fn ot_r2() -> u16 {
    FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_LAST_SEGMENT | FLAG_REVERSE
}
fn ca_total(env: &TestEnv) -> u64 {
    genome_stats(&env.stats)["CA_total"].parse().unwrap()
}

#[test]
fn long_insert_pair_trims_only_outer_ends() {
    // R1 covers ref[0,10), R2 covers ref[20,30): a long insert, mates far apart.
    // The two fragment termini are R1's 5' (ref[0,3), → C@0,2) and R2's 5'
    // (reverse, the high end ref[27,30), → C@28). The interior read ends — R1's
    // 3' (C@8) and R2's 3' (C@20) — are NOT termini and stay counted.
    let refseq = "CA".repeat(15); // 30 bp, C @ 0,2,..,28
    let reference = RefBuilder::new().contig("chr1", &refseq);
    let sam = SamBuilder::new()
        .sq("chr1", refseq.len())
        .record("p", ot_r1(), "chr1", 1, "10M", "CACACACACA", &q40(10)) // ref[0,10)
        .record("p", ot_r2(), "chr1", 21, "10M", "CACACACACA", &q40(10)); // ref[20,30)

    // Baseline: no trim → all 10 monitored C's counted.
    let env0 = TestEnv::new();
    let s0 = env0.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env0, &["--stats", &s0]);
    assert_eq!(ca_total(&env0), 10);

    // Trim 3: removes the two outer termini (C@0,2 and C@28) → 7 remain.
    let env = TestEnv::new();
    let s = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &s]);
    assert_eq!(ca_total(&env), 7, "only the two outer fragment ends are trimmed");
}

#[test]
fn overlapping_pair_trims_interior_terminus_via_genomic_skip() {
    // Fully overlapping mates over ref[0,10). Overlap dedup assigns the span to
    // R1, so R2 contributes nothing. R1's own 5' (C@0,2) is trimmed by its read
    // window; the RIGHT terminus (R2's 5' = ref[7,10)) falls in R1's interior and
    // is trimmed via the genomic mate-terminus skip — exactly the base (C@8) that
    // overlap dedup alone would have left counted. Leaves C@4,6 → 2.
    let refseq = "CACACACACA"; // C @ 0,2,4,6,8
    let reference = RefBuilder::new().contig("chr1", refseq);
    let sam = SamBuilder::new()
        .sq("chr1", refseq.len())
        .record("p", ot_r1(), "chr1", 1, "10M", "CACACACACA", &q40(10))
        .record("p", ot_r2(), "chr1", 1, "10M", "CACACACACA", &q40(10));

    // Baseline: overlap dedup → R1's 5 C's only.
    let env0 = TestEnv::new();
    let s0 = env0.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env0, &["--stats", &s0]);
    assert_eq!(ca_total(&env0), 5);

    let env = TestEnv::new();
    let s = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &s]);
    assert_eq!(ca_total(&env), 2, "both genomic termini trimmed; interior C@4,6 remain");
}

#[test]
fn orphan_read_falls_back_to_both_ends() {
    // A paired R1 whose mate is absent/unmapped: with only one mapped read the
    // far template end can't be located, so both of its ends are trimmed (like
    // single-end). keep[3,7) → C@4,6 → 2.
    let refseq = "CACACACACA";
    let reference = RefBuilder::new().contig("chr1", refseq);
    let sam = SamBuilder::new().sq("chr1", refseq.len()).record(
        "p",
        FLAG_PAIRED | FLAG_PROPER_PAIR | FLAG_FIRST_SEGMENT,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        &q40(10),
    );
    let env = TestEnv::new();
    let s = env.stats.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &s]);
    assert_eq!(ca_total(&env), 2, "orphan trims both ends");
}
