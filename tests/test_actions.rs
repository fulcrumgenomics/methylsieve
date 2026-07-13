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
