//! Actions on unconverted templates: QC-fail by default, `--no-qc-fail`,
//! `--remove-unconverted`, and a custom `--tag`.

mod helpers;
use helpers::*;

const REF: &str = "CACACACACA"; // 5 CpA → fully-C read is unconverted

fn one_unconverted_read() -> SamBuilder {
    SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "10M", "CACACACACA", &q40(10))
}

#[test]
fn default_tags_and_qc_fails() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_unconverted_read(), &reference, &env, &[]);
    assert_eq!(tag_string(&recs[0], *b"XX").as_deref(), Some("UC"));
    assert!(u16::from(recs[0].flags()) & FLAG_QC_FAIL != 0, "0x200 set by default");
}

#[test]
fn no_qc_fail_keeps_tag_without_flag() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_unconverted_read(), &reference, &env, &["--no-qc-fail"]);
    assert!(has_tag(&recs[0], *b"XX"), "tag still set");
    assert_eq!(u16::from(recs[0].flags()) & FLAG_QC_FAIL, 0, "0x200 must be off");
}

#[test]
fn remove_unconverted_drops_records() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_unconverted_read(), &reference, &env, &["--remove-unconverted"]);
    assert!(recs.is_empty(), "unconverted template should be dropped entirely");
}

#[test]
fn custom_tag_is_honored() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&one_unconverted_read(), &reference, &env, &["--tag", "YY:Z:FOO"]);
    assert_eq!(tag_string(&recs[0], *b"YY").as_deref(), Some("FOO"));
    assert!(!has_tag(&recs[0], *b"XX"), "default tag not present when overridden");
}

#[test]
fn conversion_tag_replaces_a_four_byte_foreign_tag_of_the_same_name() {
    // `--tag` takes the same replacement path as `--count-tag`, so it has the same
    // exposure: a same-named field of another type whose size matches the new
    // value's must come back as a `Z` string, not as the old type with new bytes
    // written over it. `XX:i:100000` is four bytes; `UC` plus its NUL is three.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        &q40(10),
        &["XX:i:100000"],
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(physical_tag_counts(&env.output, *b"XX"), vec![1], "exactly one XX field");
    assert_eq!(tag_string(&recs[0], *b"XX").as_deref(), Some("UC"));
}

#[test]
fn conversion_tag_from_a_previous_run_is_not_duplicated() {
    // Re-running over already-marked output must leave one XX field carrying the
    // current value, not a second copy beside the inherited one.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        &q40(10),
        &["XX:Z:UC"],
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(physical_tag_counts(&env.output, *b"XX"), vec![1], "no duplicate XX");
    assert_eq!(tag_string(&recs[0], *b"XX").as_deref(), Some("UC"));
}
