//! The three `--ref-encoding` paths must produce identical results on a pure
//! ACGT reference (where 2-bit folds nothing). This guards the default
//! (`twobit`) against drift from `bytes`/`nibble` and locks in cross-encoding
//! equivalence at the binary level.

mod helpers;
use helpers::*;

// Mixed-context, pure-ACGT reference (no N): every encoding must agree.
const REF: &str = "CACGCACGCATTACGTCACA";

#[test]
fn all_encodings_agree_on_acgt_reference() {
    let reference = RefBuilder::new().contig("chr1", REF);
    // A forward read and a reverse read so both monitored strands are exercised.
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("fwd", 0, "chr1", 1, "20M", REF, &q40(20))
        .record("rev", FLAG_REVERSE, "chr1", 1, "20M", REF, &q40(20));

    let mut rows = Vec::new();
    for enc in ["bytes", "nibble", "twobit"] {
        let env = TestEnv::new();
        let stats = env.metrics_prefix.to_str().unwrap().to_string();
        run_ok(&sam, &reference, &env, &["--ref-encoding", enc, "--metrics-prefix", &stats]);
        rows.push((enc, genome_stats(&env.stats)));
    }

    // Every per-context and summary field must match the `bytes` reference.
    let (_, base) = &rows[0];
    for (enc, row) in &rows[1..] {
        assert_eq!(row, base, "{enc} stats differ from bytes on an ACGT reference");
    }
    // Sanity: the input actually produced monitored sites.
    let ca: u64 = base["CA_total"].parse().unwrap();
    let cg: u64 = base["CG_total"].parse().unwrap();
    assert!(ca + cg > 0, "expected some monitored sites");
}

#[test]
fn default_encoding_is_twobit_equivalent_on_acgt() {
    // Running with no --ref-encoding (default twobit) must match an explicit
    // bytes run on a pure-ACGT reference.
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam =
        SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "20M", REF, &q40(20));

    let env_default = TestEnv::new();
    let sd = env_default.metrics_prefix.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env_default, &["--metrics-prefix", &sd]);

    let env_bytes = TestEnv::new();
    let sb = env_bytes.metrics_prefix.to_str().unwrap().to_string();
    run_ok(&sam, &reference, &env_bytes, &["--ref-encoding", "bytes", "--metrics-prefix", &sb]);

    assert_eq!(genome_stats(&env_default.stats), genome_stats(&env_bytes.stats));
}
