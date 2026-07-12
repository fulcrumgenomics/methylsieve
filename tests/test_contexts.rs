//! `--contexts` selects which contexts count toward the threshold, flipping
//! the decision for the same input.

mod helpers;
use helpers::*;

// C at 0,2,4,6,8: contexts CpA(0), CpG(2), CpA(4), CpG(6), CpA(8) → 3 CpA, 2 CpG.
const REF: &str = "CACGCACGCA";

fn one_read() -> SamBuilder {
    SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "10M", REF, &q40(10))
}

#[test]
fn default_cph_counts_three_cpa_and_tags() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_read(), &reference, &env, &[]);
    assert!(has_tag(&recs[0], *b"XX"), "3 CpA ≥ 3 under default CpH → tagged");
}

#[test]
fn cpg_only_counts_two_and_does_not_tag() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_read(), &reference, &env, &["--contexts", "CpG"]);
    assert!(!has_tag(&recs[0], *b"XX"), "2 CpG < 3 → not tagged");
}

#[test]
fn cpa_plus_cpg_counts_five_and_tags() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_read(), &reference, &env, &["--contexts", "CpA,CpG"]);
    assert!(has_tag(&recs[0], *b"XX"), "3 CpA + 2 CpG = 5 ≥ 3 → tagged");
}
