//! The `@PG` line is appended and chains correctly across re-runs.

mod helpers;
use helpers::*;

use std::process::Command;

const REF: &str = "CACACACACA";

#[test]
fn pg_line_is_appended_and_chains_on_rerun() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam =
        SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "10M", REF, &q40(10));

    // First run: SAM in, BAM out.
    let r1 = run_methylsieve(&sam, &reference, &env, &[]);
    assert!(r1.status_ok, "first run failed: {}", r1.stderr);
    let pg1 = read_pg_lines(&env.output);
    assert_eq!(pg1.iter().filter(|l| l.contains("PN:methylsieve")).count(), 1);

    // Second run: feed the first BAM back through methylsieve.
    let second = env._tmp.path().join("out2.bam");
    let out = Command::new(methylsieve_binary())
        .arg("-i")
        .arg(&env.output)
        .arg("-o")
        .arg(&second)
        .arg("-r")
        .arg(&env.reference)
        .output()
        .expect("second run");
    assert!(out.status.success(), "second run failed: {}", String::from_utf8_lossy(&out.stderr));

    let pg2 = read_pg_lines(&second);
    assert_eq!(
        pg2.iter().filter(|l| l.contains("PN:methylsieve")).count(),
        2,
        "two methylsieve @PG records after two runs: {pg2:?}"
    );
    assert!(pg2.iter().any(|l| l.contains("PP:")), "second @PG must chain via PP:");
}
