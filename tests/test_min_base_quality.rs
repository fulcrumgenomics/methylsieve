//! Read bases below `--min-base-quality` are not tallied.

mod helpers;
use helpers::*;

const REF: &str = "CACACACACA"; // 5 CpA C's at 0,2,4,6,8

/// Quals: Phred 2 (`#`) at the C positions, Phred 40 (`I`) elsewhere.
const LOW_C_QUALS: &str = "#I#I#I#I#I";

#[test]
fn low_quality_cytosines_are_skipped() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        LOW_C_QUALS,
    );
    // Default --min-base-quality 10 skips the BQ-2 C's → 0 unconverted.
    let recs = run_ok(&sam, &reference, &env, &[]);
    assert!(!has_tag(&recs[0], *b"XX"), "BQ-2 cytosines must not count");
}

#[test]
fn min_base_quality_zero_counts_them() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        LOW_C_QUALS,
    );
    let recs = run_ok(&sam, &reference, &env, &["--min-base-quality", "0"]);
    assert!(has_tag(&recs[0], *b"XX"), "with -q0 the 5 unconverted C's tag the read");
}
