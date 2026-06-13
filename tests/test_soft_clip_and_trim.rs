//! Soft-clip handling and `--ignore-template-ends` trimming for **single-end**
//! reads, verified through the per-context tallies in the `--stats` TSV. A lone
//! read's far template terminus can't be located, so both of its ends are
//! trimmed. Paired-end genomic-terminus behavior is covered in
//! `test_template_ends.rs`.

mod helpers;
use helpers::*;

const TOP: &str = "CACACACACA"; // 5 CpA C's at 0,2,4,6,8
const BOT: &str = "TGTGTGTGTG"; // 5 CpA-equivalent G's at 1,3,5,7,9

fn ca_unconv(env: &TestEnv) -> u64 {
    genome_stats(&env.stats)["CA_unconv"].parse().unwrap()
}
fn ca_total(env: &TestEnv) -> u64 {
    genome_stats(&env.stats)["CA_total"].parse().unwrap()
}

#[test]
fn soft_clipped_bases_are_not_tallied() {
    // 3S7M: the 3 clipped bases never align; the 7M covers ref 0..6, hitting
    // C@0,2,4,6 (all read as C) → 4 unconverted, 4 total. No trimming applied.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", TOP);
    let stats = env.stats.to_str().unwrap().to_string();
    let sam = SamBuilder::new().sq("chr1", TOP.len()).record(
        "r",
        0,
        "chr1",
        1,
        "3S7M",
        "GGGCACACAC",
        &q40(10),
    );
    run_ok(&sam, &reference, &env, &["--stats", &stats]);
    assert_eq!(ca_total(&env), 4);
    assert_eq!(ca_unconv(&env), 4);
}

#[test]
fn single_end_trims_both_ends_forward() {
    // A single-end read trims BOTH ends: --ignore-template-ends 3 keeps stored
    // positions [3,7). Drops C@0,C@2 (5') and C@8 (3'); leaves C@4,C@6 → 2.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", TOP);
    let stats = env.stats.to_str().unwrap().to_string();
    let sam = SamBuilder::new().sq("chr1", TOP.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        &q40(10),
    );
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &stats]);
    assert_eq!(ca_total(&env), 2);
    assert_eq!(ca_unconv(&env), 2);
}

#[test]
fn single_end_trims_both_ends_reverse() {
    // A reverse single-end read on a ref of G's monitors G@1,3,5,7,9. Both ends
    // are trimmed over *stored* positions, keeping [3,7): drops G@1 (low) and
    // G@7,G@9 (high); leaves G@3,G@5 → 2 unconverted.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", BOT);
    let stats = env.stats.to_str().unwrap().to_string();
    let sam = SamBuilder::new().sq("chr1", BOT.len()).record(
        "r",
        FLAG_REVERSE,
        "chr1",
        1,
        "10M",
        "TGTGTGTGTG",
        &q40(10),
    );
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &stats]);
    assert_eq!(ca_total(&env), 2);
    assert_eq!(ca_unconv(&env), 2);
}

#[test]
fn soft_clip_counts_toward_trim_budget() {
    // 3S7M with --ignore-template-ends 3: the 3 leading soft-clips absorb the 5'
    // budget (no aligned base dropped there), while the 3' trim drops stored
    // 7,8,9 (= C@ref4,ref6). Aligned C's at ref0,2,4,6 → keep ref0,ref2 → 2.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", TOP);
    let stats = env.stats.to_str().unwrap().to_string();
    let sam = SamBuilder::new().sq("chr1", TOP.len()).record(
        "r",
        0,
        "chr1",
        1,
        "3S7M",
        "GGGCACACAC",
        &q40(10),
    );
    run_ok(&sam, &reference, &env, &["--ignore-template-ends", "3", "--stats", &stats]);
    assert_eq!(ca_unconv(&env), 2, "5' clip absorbs the 5' budget; 3' trim still drops 2");
}
