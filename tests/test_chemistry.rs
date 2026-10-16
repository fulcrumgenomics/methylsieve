//! `--chemistry`: chemistries that convert methylated cytosines (taps,
//! illumina-5base) make no unconverted call and read M-bias as methylation, and
//! the chemistry check stops a run whose data plainly belongs to the other kind.
//!
//! Most fixtures tile the reference with `CACG`, so a forward read monitors one
//! CpA and one CpG cytosine per 4 bp; each read converts a chosen number of each.
//! The chemistry check only stops a run past 100k CpH and 10k CpG sites, so the
//! fixtures that exercise it use enough reads to clear both floors.

mod helpers;
use helpers::*;

/// Number of `CACG` units in the reference (and in every read).
const UNITS: usize = 25;

/// Reads per fixture that clear the chemistry check's site floors: each read
/// carries `UNITS` CpA and `UNITS` CpG sites, so this gives 102.5k of each.
const READS_PAST_FLOORS: usize = 4100;

/// The `CACG`-tiled reference.
fn cacg_ref() -> String {
    "CACG".repeat(UNITS)
}

/// A forward read over the reference whose first `cph_converted` CpA cytosines
/// and first `cpg_converted` CpG cytosines read as T; the rest read as C.
fn cacg_read(cph_converted: usize, cpg_converted: usize) -> String {
    (0..UNITS)
        .map(|u| {
            let cpa = if u < cph_converted { 'T' } else { 'C' };
            let cpg = if u < cpg_converted { 'T' } else { 'C' };
            format!("{cpa}A{cpg}G")
        })
        .collect()
}

/// `n` single-end forward templates, all with the same conversion pattern.
fn cacg_sam(n: usize, cph_converted: usize, cpg_converted: usize) -> SamBuilder {
    let len = UNITS * 4;
    let mut b = SamBuilder::new().sq("chr1", len);
    let seq = cacg_read(cph_converted, cpg_converted);
    let qual = q40(len);
    for i in 0..n {
        b = b.record(&format!("r{i}"), 0, "chr1", 1, &format!("{len}M"), &seq, &qual);
    }
    b
}

fn cacg_reference() -> RefBuilder {
    RefBuilder::new().contig("chr1", &cacg_ref())
}

/// TAPS-like reads: no CpH converted, about half of CpG converted.
fn taps_like(n: usize) -> SamBuilder {
    cacg_sam(n, 0, 12)
}

/// EM-seq-like reads: every CpH converted, a quarter of CpG converted.
fn emseq_like(n: usize) -> SamBuilder {
    cacg_sam(n, UNITS, 6)
}

/// Reads in which nothing converted.
fn unconverted(n: usize) -> SamBuilder {
    cacg_sam(n, 0, 0)
}

/// CpG units at the start of [`taps_like_with_converted_start`] reads that all
/// convert.
const CONVERTED_START_UNITS: usize = 5;

/// TAPS-like reads that, read as em-seq, dip at the 5' end: the first
/// [`CONVERTED_START_UNITS`] CpGs all convert, then read `i` converts every other
/// CpG, alternating which by read, so the rest of the curve sits at one half. No
/// CpH converts. Masking the dip drops a fifth of each read's CpA sites.
fn taps_like_with_converted_start(n: usize) -> SamBuilder {
    let len = UNITS * 4;
    let mut b = SamBuilder::new().sq("chr1", len);
    let qual = q40(len);
    for i in 0..n {
        let seq: String = (0..UNITS)
            .map(|u| {
                let converted = u < CONVERTED_START_UNITS || (u + i) % 2 == 0;
                format!("CA{}G", if converted { 'T' } else { 'C' })
            })
            .collect();
        b = b.record(&format!("r{i}"), 0, "chr1", 1, &format!("{len}M"), &seq, &qual);
    }
    b
}

// ── No unconverted call ─────────────────────────────────────────────────────

#[test]
fn taps_never_tags_or_qc_fails_a_fully_retained_template() {
    let env = TestEnv::new();
    let recs = run_ok(&unconverted(10), &cacg_reference(), &env, &["--chemistry", "taps"]);
    assert_eq!(recs.len(), 10);
    for rec in &recs {
        assert_eq!(u16::from(rec.flags()) & FLAG_QC_FAIL, 0);
        assert!(!has_tag(rec, *b"XX"), "no conversion tag under taps");
        assert!(!has_tag(rec, *b"ch"), "no count tag under taps");
    }
}

#[test]
fn emseq_flags_the_same_fully_retained_templates() {
    // The contrast for the test above: the same reads are unconverted under em-seq.
    let env = TestEnv::new();
    let recs = run_ok(&unconverted(10), &cacg_reference(), &env, &[]);
    assert!(recs.iter().all(|r| u16::from(r.flags()) & FLAG_QC_FAIL != 0));
}

#[test]
fn taps_keeps_all_templates_despite_remove_unconverted() {
    let env = TestEnv::new();
    let recs = run_ok(
        &unconverted(10),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps", "--remove-unconverted"],
    );
    assert_eq!(recs.len(), 10);
}

#[test]
fn taps_warns_listing_explicit_decision_flags() {
    let env = TestEnv::new();
    let r = run_methylsieve(
        &unconverted(10),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps", "--remove-unconverted", "--min-sites", "10"],
    );
    assert!(r.status_ok, "{}", r.stderr);
    assert!(
        r.stderr.contains("--min-sites, --remove-unconverted have no effect with --chemistry taps"),
        "{}",
        r.stderr
    );
}

#[test]
fn taps_does_not_warn_about_defaulted_decision_flags() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    let r = run_methylsieve(
        &unconverted(10),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps", "--metrics-prefix", &prefix],
    );
    assert!(r.status_ok, "{}", r.stderr);
    assert!(!r.stderr.contains("no effect"), "{}", r.stderr);
    assert!(!r.stderr.contains("pass through unchanged"), "{}", r.stderr);
}

#[test]
fn taps_without_mask_or_metrics_passes_input_through_with_warning() {
    let env = TestEnv::new();
    let r = run_methylsieve(&unconverted(10), &cacg_reference(), &env, &["--chemistry", "taps"]);
    assert!(r.status_ok, "{}", r.stderr);
    assert!(r.stderr.contains("records pass through unchanged"), "{}", r.stderr);
    assert_eq!(read_records(&r.output).len(), 10);
}

#[test]
fn taps_without_controls_writes_no_conversion_matrix() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    run_ok(
        &taps_like(10),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps", "--metrics-prefix", &prefix],
    );
    assert!(env.stats.exists(), "summary.tsv is still written");
    let matrix = std::path::PathBuf::from(format!("{prefix}.conversion-matrix.tsv"));
    assert!(!matrix.exists(), "no conversion matrix without an unconverted call");
}

#[test]
fn taps_summary_records_chemistry_and_methylation() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    run_ok(
        &taps_like(10),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps", "--metrics-prefix", &prefix],
    );
    let row = genome_stats(&env.stats);
    assert_eq!(row["chemistry"], "taps");
    assert_eq!(row["n_unconverted"], "");
    // 12 of 25 CpGs converted, which TAPS reads as methylated.
    assert_eq!(row["CpG_meth_rate"], "0.480000");
}

// ── Chemistries that share a readout ────────────────────────────────────────

#[test]
fn illumina_5base_output_matches_taps() {
    let (taps_env, i5b_env) = (TestEnv::new(), TestEnv::new());
    let taps = run_ok(&unconverted(10), &cacg_reference(), &taps_env, &["--chemistry", "taps"]);
    let i5b =
        run_ok(&unconverted(10), &cacg_reference(), &i5b_env, &["--chemistry", "illumina-5base"]);
    assert_eq!(taps, i5b);
}

#[test]
fn bisulfite_output_matches_emseq() {
    let (emseq_env, bis_env) = (TestEnv::new(), TestEnv::new());
    let emseq = run_ok(&unconverted(10), &cacg_reference(), &emseq_env, &[]);
    let bisulfite =
        run_ok(&unconverted(10), &cacg_reference(), &bis_env, &["--chemistry", "bisulfite"]);
    assert_eq!(emseq, bisulfite);
}

// ── M-bias in methylation space ─────────────────────────────────────────────

/// CpG-dense reference: a CpG at every even position.
const CPG_REF: &str = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";

/// `n` forward reads over [`CPG_REF`] whose first 8 cycles keep their CpG
/// cytosines as C and convert the rest to T. Under TAPS that is a methylation
/// dip at the read start, the shape end-repair fill-in leaves.
fn retained_start_sam(n: usize) -> SamBuilder {
    let seq: String = (0..CPG_REF.len())
        .map(|i| {
            if i % 2 == 1 {
                'G'
            } else if i < 8 {
                'C'
            } else {
                'T'
            }
        })
        .collect();
    let mut b = SamBuilder::new().sq("chr1", CPG_REF.len());
    let qual = q40(CPG_REF.len());
    for i in 0..n {
        b = b.record(&format!("r{i}"), 0, "chr1", 1, "40M", &seq, &qual);
    }
    b
}

#[test]
fn taps_masks_a_retained_read_start() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", CPG_REF);
    let recs = run_ok(
        &retained_start_sam(150),
        &reference,
        &env,
        &["--chemistry", "taps", "--mbias-mask"],
    );
    let masked = leading_quality_run(&recs[0], 2);
    assert!((5..=15).contains(&masked), "expected ~8 leading Q2 cycles, got {masked}");
}

#[test]
fn emseq_does_not_mask_the_same_retained_read_start() {
    // Read as em-seq, the retained start is *high* methylation, not a dip.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", CPG_REF);
    let recs = run_ok(&retained_start_sam(150), &reference, &env, &["--mbias-mask"]);
    assert_eq!(leading_quality_run(&recs[0], 2), 0);
}

// ── Chemistry check ─────────────────────────────────────────────────────────

#[test]
fn emseq_chemistry_on_taps_like_data_errors() {
    let env = TestEnv::new();
    let r = run_methylsieve(&taps_like(READS_PAST_FLOORS), &cacg_reference(), &env, &[]);
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
}

#[test]
fn taps_chemistry_on_emseq_like_data_errors() {
    let env = TestEnv::new();
    let r = run_methylsieve(
        &emseq_like(READS_PAST_FLOORS),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps"],
    );
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry em-seq"), "{}", r.stderr);
}

#[test]
fn mismatch_is_caught_when_input_ends_while_learning_the_mask() {
    let env = TestEnv::new();
    let r =
        run_methylsieve(&taps_like(READS_PAST_FLOORS), &cacg_reference(), &env, &["--mbias-mask"]);
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
}

#[test]
fn mismatch_is_caught_when_the_mask_buffer_fills() {
    // The full buffer alone clears the floors, so the check runs at the freeze,
    // on exactly the buffered templates' sites.
    let env = TestEnv::new();
    let buffer = READS_PAST_FLOORS.to_string();
    let r = run_methylsieve(
        &taps_like(2 * READS_PAST_FLOORS),
        &cacg_reference(),
        &env,
        &["--mbias-mask", "--mbias-buffer-templates", &buffer],
    );
    assert!(!r.status_ok);
    let buffered_sites = READS_PAST_FLOORS * UNITS;
    assert!(r.stderr.contains(&format!("CpH conversion 0.0000 over {buffered_sites} sites")));
}

#[test]
fn mismatch_is_caught_while_streaming_past_a_small_mask_buffer() {
    let env = TestEnv::new();
    let r = run_methylsieve(
        &taps_like(READS_PAST_FLOORS),
        &cacg_reference(),
        &env,
        &["--mbias-mask", "--mbias-buffer-templates", "100"],
    );
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
}

#[test]
fn masked_sites_still_count_toward_the_chemistry_check() {
    // The learned mask removes enough CpA sites from the tally to leave it below
    // 100k, but the check counts every site seen before masking.
    let env = TestEnv::new();
    let r = run_methylsieve(
        &taps_like_with_converted_start(READS_PAST_FLOORS),
        &cacg_reference(),
        &env,
        &["--mbias-mask", "--mbias-buffer-templates", "100"],
    );
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
}

#[test]
fn hypomethylated_taps_data_declared_emseq_errors() {
    // 3 of 25 CpGs converted (12%), no CpH: still plainly TAPS-like.
    let env = TestEnv::new();
    let r = run_methylsieve(&cacg_sam(READS_PAST_FLOORS, 0, 3), &cacg_reference(), &env, &[]);
    assert!(!r.status_ok);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
}

#[test]
fn entirely_unconverted_data_declared_emseq_passes_check() {
    // A failed em-seq library: the check stays out of the way and every template
    // is flagged unconverted, which is the point of em-seq mode.
    let env = TestEnv::new();
    let recs = run_ok(&unconverted(READS_PAST_FLOORS), &cacg_reference(), &env, &[]);
    assert!(recs.iter().all(|r| u16::from(r.flags()) & FLAG_QC_FAIL != 0));
}

#[test]
fn entirely_unconverted_data_declared_taps_passes_check() {
    let env = TestEnv::new();
    let recs =
        run_ok(&unconverted(READS_PAST_FLOORS), &cacg_reference(), &env, &["--chemistry", "taps"]);
    assert_eq!(recs.len(), READS_PAST_FLOORS);
}

#[test]
fn partly_failed_emseq_library_passes_check() {
    // Half the CpH converted: neither pattern, so no verdict either way.
    let env = TestEnv::new();
    let r = run_methylsieve(&cacg_sam(READS_PAST_FLOORS, 12, 3), &cacg_reference(), &env, &[]);
    assert!(r.status_ok, "{}", r.stderr);
}

#[test]
fn ambiguous_data_logs_class_without_erroring() {
    let env = TestEnv::new();
    let r = run_methylsieve(
        &cacg_sam(READS_PAST_FLOORS, 12, 3),
        &cacg_reference(),
        &env,
        &["--chemistry", "taps"],
    );
    assert!(r.status_ok, "{}", r.stderr);
    assert!(r.stderr.contains("ambiguous, so not verified"), "{}", r.stderr);
}

#[test]
fn chemistry_check_below_site_floors_warns_instead_of_erroring() {
    let env = TestEnv::new();
    let r = run_methylsieve(&taps_like(10), &cacg_reference(), &env, &[]);
    assert!(r.status_ok, "{}", r.stderr);
    assert!(r.stderr.contains("Re-run with --chemistry taps"), "{}", r.stderr);
    assert!(r.stderr.contains("Too few sites to be certain"), "{}", r.stderr);
}

#[test]
fn spike_in_controls_are_excluded_from_the_chemistry_check() {
    // A TAPS-like control twice the size of an EM-seq-like genome: counted with
    // the genome, it would make the data look ambiguous.
    let env = TestEnv::new();
    let r = run_methylsieve(
        &genome_and_control_sam(READS_PAST_FLOORS, (UNITS, 6), 2 * READS_PAST_FLOORS, (0, 12)),
        &genome_and_control_reference(),
        &env,
        &["--control-contig", "ctrl"],
    );
    assert!(r.status_ok, "{}", r.stderr);
    assert!(r.stderr.contains("EM-seq-like, consistent"), "{}", r.stderr);
}

// ── Run statistics ──────────────────────────────────────────────────────────

#[test]
fn taps_run_stats_report_genome_methylation() {
    let env = TestEnv::new();
    let r = run_methylsieve(&taps_like(10), &cacg_reference(), &env, &["--chemistry", "taps"]);
    assert!(r.status_ok, "{}", r.stderr);
    assert!(
        r.stderr.contains(
            "Genome CpG apparent methylation 0.48000 over 250 sites; CpH conversion 0.00000 \
             over 250 sites."
        ),
        "{}",
        r.stderr
    );
}

#[test]
fn taps_run_stats_report_na_without_genome_sites() {
    let env = TestEnv::new();
    let r = run_methylsieve(
        &genome_and_control_sam(0, (0, 0), 3, (0, 20)),
        &genome_and_control_reference(),
        &env,
        &["--chemistry", "taps", "--control-contig", "ctrl"],
    );
    assert!(r.status_ok, "{}", r.stderr);
    assert!(
        r.stderr.contains(
            "Genome CpG apparent methylation n/a over 0 sites; CpH conversion n/a over 0 sites."
        ),
        "{}",
        r.stderr
    );
}

// ── Control conversion matrices ─────────────────────────────────────────────

/// A SAM over a `chr1` genome contig and a `ctrl` control contig, both `CACG`
/// tiled: `n_genome` genome templates and `n_ctrl` control templates, each with
/// the given numbers of CpA and CpG cytosines converted.
fn genome_and_control_sam(
    n_genome: usize,
    genome: (usize, usize),
    n_ctrl: usize,
    ctrl: (usize, usize),
) -> SamBuilder {
    let len = UNITS * 4;
    let qual = q40(len);
    let cigar = format!("{len}M");
    let mut b = SamBuilder::new().sq("chr1", len).sq("ctrl", len);
    for i in 0..n_genome {
        b = b.record(&format!("g{i}"), 0, "chr1", 1, &cigar, &cacg_read(genome.0, genome.1), &qual);
    }
    for i in 0..n_ctrl {
        b = b.record(&format!("c{i}"), 0, "ctrl", 1, &cigar, &cacg_read(ctrl.0, ctrl.1), &qual);
    }
    b
}

fn genome_and_control_reference() -> RefBuilder {
    RefBuilder::new().contig("chr1", &cacg_ref()).contig("ctrl", &cacg_ref())
}

/// The rows of `PREFIX.conversion-matrix.tsv`.
fn matrix_rows(env: &TestEnv) -> Vec<std::collections::HashMap<String, String>> {
    read_stats_rows(&std::path::PathBuf::from(format!(
        "{}.conversion-matrix.tsv",
        env.metrics_prefix_arg()
    )))
}

#[test]
fn taps_control_matrix_counts_cpg_without_verdicts() {
    // A CpG-methylated control: 20 of 25 CpGs converted, no CpH.
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    run_ok(
        &genome_and_control_sam(10, (0, 12), 3, (0, 20)),
        &genome_and_control_reference(),
        &env,
        &["--chemistry", "taps", "--control-contig", "ctrl", "--metrics-prefix", &prefix],
    );
    let rows = matrix_rows(&env);
    assert_eq!(rows.len(), 1, "only the control has a matrix under taps: {rows:?}");
    let row = &rows[0];
    assert_eq!(row["scope"], "ctrl");
    assert_eq!(row["context"], "CpG");
    assert_eq!(row["checked_sites"], "25");
    assert_eq!(row["converted_sites"], "20");
    assert_eq!(row["n_templates"], "3");
    assert_eq!(row["decision"], "");
    assert!(std::path::PathBuf::from(format!("{prefix}.conversion-matrix.pdf")).exists());
}

#[test]
fn taps_summary_leaves_control_unconverted_blank() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    run_ok(
        &genome_and_control_sam(10, (0, 12), 3, (0, 0)),
        &genome_and_control_reference(),
        &env,
        &["--chemistry", "taps", "--control-contig", "ctrl", "--metrics-prefix", &prefix],
    );
    let rows = read_stats_rows(&env.stats);
    let ctrl = rows.iter().find(|r| r["scope"] == "ctrl").unwrap();
    assert_eq!(ctrl["n_unconverted"], "");
}

#[test]
fn emseq_control_matrix_follows_genome_with_verdicts() {
    // Converted genome templates; control templates that failed to convert.
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    let recs = run_ok(
        &genome_and_control_sam(10, (UNITS, 6), 3, (0, 0)),
        &genome_and_control_reference(),
        &env,
        &["--control-contig", "ctrl", "--metrics-prefix", &prefix],
    );
    let rows = matrix_rows(&env);
    let scopes: Vec<&str> = rows.iter().map(|r| r["scope"].as_str()).collect();
    assert_eq!(scopes, ["genome", "ctrl"]);
    let ctrl = &rows[1];
    assert_eq!(ctrl["context"], "CpH");
    assert_eq!(ctrl["converted_sites"], "0");
    assert_eq!(ctrl["decision"], "unconverted");
    // Counted against the boundary, but control records are never flagged.
    assert!(recs.iter().all(|r| u16::from(r.flags()) & FLAG_QC_FAIL == 0));
}

#[test]
fn emseq_summary_counts_unconverted_control_templates() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    run_ok(
        &genome_and_control_sam(10, (UNITS, 6), 3, (0, 0)),
        &genome_and_control_reference(),
        &env,
        &["--control-contig", "ctrl", "--metrics-prefix", &prefix],
    );
    let rows = read_stats_rows(&env.stats);
    let ctrl = rows.iter().find(|r| r["scope"] == "ctrl").unwrap();
    assert_eq!(ctrl["n_unconverted"], "3");
    assert_eq!(ctrl["frac_unconverted"], "1.000000");
}

/// Pages in a PDF: page objects are written as `/Type/Page`, the page tree as
/// `/Type/Pages`.
fn pdf_page_count(path: &std::path::Path) -> usize {
    let bytes = std::fs::read(path).expect("PDF written");
    let marker = b"/Type/Page";
    bytes
        .windows(marker.len() + 1)
        .filter(|w| w.starts_with(marker) && w[marker.len()] != b's')
        .count()
}

#[test]
fn emseq_conversion_matrix_pdf_has_a_page_per_scope() {
    let env = TestEnv::new();
    let prefix = env.metrics_prefix_arg();
    let len = UNITS * 4;
    let (qual, cigar) = (q40(len), format!("{len}M"));
    let mut sam = SamBuilder::new().sq("chr1", len).sq("lambda", len).sq("puc19", len);
    for (contig, seq) in
        [("chr1", cacg_read(UNITS, 6)), ("lambda", cacg_read(UNITS, 0)), ("puc19", cacg_read(0, 0))]
    {
        for i in 0..3 {
            sam = sam.record(&format!("{contig}{i}"), 0, contig, 1, &cigar, &seq, &qual);
        }
    }
    let reference = RefBuilder::new()
        .contig("chr1", &cacg_ref())
        .contig("lambda", &cacg_ref())
        .contig("puc19", &cacg_ref());
    let args =
        ["--control-contig", "lambda", "--control-contig", "puc19", "--metrics-prefix", &prefix];
    run_ok(&sam, &reference, &env, &args);
    let pdf = std::path::PathBuf::from(format!("{prefix}.conversion-matrix.pdf"));
    assert_eq!(pdf_page_count(&pdf), 3);
}
