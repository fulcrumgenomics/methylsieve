//! Metric TSVs written under `--metrics-prefix PREFIX`.
//!
//! Each file is one serde row struct serialized with [`fgoxide::io::DelimFile`]
//! (the same pattern as riker): column names and order come from the struct's
//! fields, and rates are rendered as fixed 6-decimal fractions (never
//! percentages) via the `ser_*` helpers.
//!
//! - `PREFIX.summary.tsv` — one per-context conversion row per scope (genome
//!   first, then each `--control-contig`). Every row is decision-basis:
//!   overlap-deduped, end-trimmed, and including supplementary evidence — i.e.
//!   exactly the sites the unconverted call acted on. The applied per-read mask
//!   lengths (`r1_mask_5p`/`r2_mask_5p`/`se_mask_5p`/`se_mask_3p`, in sequencing
//!   cycles) ride along as run-level columns, blank when masking was not run.
//!   Per-read conversion broken out by cycle lives in `mbias.tsv`, not here.
//! - `PREFIX.conversion-matrix.tsv` — per-`(checked, converted)` template cell, for
//!   the genome (only for chemistries that make the unconverted call) and for
//!   each `--control-contig`.
//! - `PREFIX.mbias.tsv` — per-read-cycle methylation by `(read, end, context)`.
//!
//! Conversion rates are chemistry-neutral (the fraction of cytosines read as T);
//! methylation columns are read through the chemistry's
//! [`crate::chemistry::Readout`].

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow};
use fgoxide::io::DelimFile;
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;
use serde::Serialize;

use crate::METHYLSIEVE_BUILD;
use crate::chemistry::Chemistry;
use crate::mask::MaskPlan;
use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
use crate::plots::MatrixPage;
use crate::reference::Context;
use crate::sieve::{DecidedBy, PerContextCounters, ScopeStats, Stats};

/// Build the path `PREFIX.<suffix>` from a metrics prefix.
fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(".");
    s.push(suffix);
    PathBuf::from(s)
}

/// Write every metric file under `prefix`. `mbias` holds the per-cycle counts
/// accumulated during the run; `classify(unconv, checked)` replays the per-cell
/// decision verdict for the conversion matrix. `contexts` and `control_contexts`
/// label the contexts the genome's and the controls' matrices were counted over.
#[allow(clippy::too_many_arguments)] // a single wiring call site; a params struct would not aid clarity
pub(crate) fn write_all<F>(
    prefix: &Path,
    chemistry: Chemistry,
    stats: &Stats,
    mbias: &MbiasAccumulator,
    mask_plan: Option<&MaskPlan>,
    header: &Header,
    sample_override: Option<&str>,
    input_path: Option<&Path>,
    contexts: &str,
    control_contexts: &str,
    classify: F,
) -> Result<()>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    let sample = resolve_sample(header, sample_override, input_path);
    let df = DelimFile::default();
    write_tsv(
        &df,
        &with_suffix(prefix, "summary.tsv"),
        summary_rows(stats, chemistry, mask_plan, &sample),
    )?;
    let calls_unconverted = chemistry.readout().calls_unconverted();
    let pages = matrix_pages(stats, calls_unconverted, contexts, control_contexts);
    if !pages.is_empty() {
        write_tsv(
            &df,
            &with_suffix(prefix, "conversion-matrix.tsv"),
            matrix_rows(&pages, &sample, &classify),
        )?;
    }
    write_tsv(&df, &with_suffix(prefix, "mbias.tsv"), mbias_rows(mbias, &sample))?;
    // PDF plots of the same data (M-bias curves + one conversion-matrix hexbin
    // page per scope).
    crate::plots::write_mbias_pdf(
        &with_suffix(prefix, "mbias.pdf"),
        mbias,
        mask_plan,
        &sample,
        chemistry,
    )?;
    crate::plots::write_matrix_pdf(
        &with_suffix(prefix, "conversion-matrix.pdf"),
        &pages,
        &classify,
        &sample,
    )?;
    Ok(())
}

/// Write `rows` as a tab-separated file at `path`, with path context on error.
fn write_tsv<T: Serialize>(df: &DelimFile, path: &Path, rows: Vec<T>) -> Result<()> {
    df.write_tsv(path, rows.iter()).map_err(|e| anyhow!("writing {}: {e}", path.display()))
}

// ── serde rate formatting ─────────────────────────────────────────────────────

/// Serialize an `f64` as a fixed 6-decimal fraction (the metric convention).
fn ser_f64_6dp<S: serde::Serializer>(v: &f64, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_str(&format!("{v:.6}"))
}

/// Serialize `Option<f64>` as a 6-decimal fraction, or the empty string for
/// `None` — used where a rate is undefined because no sites were observed.
fn ser_opt_f64_6dp<S: serde::Serializer>(v: &Option<f64>, s: S) -> Result<S::Ok, S::Error> {
    match v {
        Some(x) => s.serialize_str(&format!("{x:.6}")),
        None => s.serialize_str(""),
    }
}

/// Conversion rate `1 - unconv/total`, or `None` when no sites were observed.
fn conv_rate(unconv: u64, total: u64) -> Option<f64> {
    (total > 0).then(|| 1.0 - unconv as f64 / total as f64)
}

/// Combined CpH `(unconv, total)` across CpA/CpC/CpT.
fn cph_counts(c: &PerContextCounters) -> (u64, u64) {
    [Context::CpA, Context::CpC, Context::CpT]
        .iter()
        .fold((0, 0), |(u, t), &ctx| (u + c.unconv_for(ctx), t + c.total_for(ctx)))
}

// ── summary.tsv ───────────────────────────────────────────────────────────────

/// One summary row — one per scope, all decision-basis. Column names/order are
/// the serde field names (CpX fields renamed to preserve the context casing).
/// The unconverted columns are blank for chemistries that make no unconverted call.
#[derive(Serialize)]
struct SummaryRow {
    sample: String,
    methylsieve_version: &'static str,
    chemistry: &'static str,
    scope: String,
    r1_mask_5p: Option<usize>,
    r2_mask_5p: Option<usize>,
    se_mask_5p: Option<usize>,
    se_mask_3p: Option<usize>,
    n_templates: u64,
    n_mapped: u64,
    n_evaluated: u64,
    n_unconverted: Option<u64>,
    #[serde(serialize_with = "ser_opt_f64_6dp")]
    frac_unconverted: Option<f64>,
    chimeric_to_control_templates: Option<u64>,
    #[serde(rename = "CpA_obs")]
    cpa_obs: u64,
    #[serde(rename = "CpA_conv_rate", serialize_with = "ser_opt_f64_6dp")]
    cpa_conv_rate: Option<f64>,
    #[serde(rename = "CpC_obs")]
    cpc_obs: u64,
    #[serde(rename = "CpC_conv_rate", serialize_with = "ser_opt_f64_6dp")]
    cpc_conv_rate: Option<f64>,
    #[serde(rename = "CpT_obs")]
    cpt_obs: u64,
    #[serde(rename = "CpT_conv_rate", serialize_with = "ser_opt_f64_6dp")]
    cpt_conv_rate: Option<f64>,
    #[serde(rename = "CpH_obs")]
    cph_obs: u64,
    #[serde(rename = "CpH_conv_rate", serialize_with = "ser_opt_f64_6dp")]
    cph_conv_rate: Option<f64>,
    #[serde(rename = "CpG_obs")]
    cpg_obs: u64,
    #[serde(rename = "CpG_conv_rate", serialize_with = "ser_opt_f64_6dp")]
    cpg_conv_rate: Option<f64>,
    // CpG methylation rate, read through the chemistry: 1 - CpG conv rate when
    // unmethylated cytosines convert, the CpG conv rate itself when methylated
    // ones do. Stated directly because it is the headline biological readout.
    #[serde(rename = "CpG_meth_rate", serialize_with = "ser_opt_f64_6dp")]
    cpg_meth_rate: Option<f64>,
}

/// The four applied mask lengths (sequencing cycles). R1/R2 are 5'-only — their
/// 3' end is the mate's domain — so only single-end reads carry a 3' mask. All
/// `None` when masking was not run; for a single-library-type file the unused
/// pair is simply blank (e.g. `se_*` for a purely paired run).
#[derive(Clone, Copy, Default)]
struct MaskCols {
    r1_5p: Option<usize>,
    r2_5p: Option<usize>,
    se_5p: Option<usize>,
    se_3p: Option<usize>,
}

impl MaskCols {
    fn from_plan(plan: Option<&MaskPlan>) -> Self {
        match plan {
            None => Self::default(),
            Some(p) => Self {
                r1_5p: Some(p.role_5p(ReadRole::R1)),
                r2_5p: Some(p.role_5p(ReadRole::R2)),
                se_5p: Some(p.role_5p(ReadRole::Se)),
                se_3p: Some(p.k_se_3p()),
            },
        }
    }
}

/// One scope row (decision basis). `chimeric` is the genome-only control
/// diagnostic (`None` for control scopes); `masks` is the run-level mask plan,
/// repeated on every row.
fn summary_row(
    sample: &str,
    chemistry: Chemistry,
    scope: &ScopeStats,
    chimeric: Option<u64>,
    masks: MaskCols,
) -> SummaryRow {
    let c = &scope.counters;
    let (cph_u, cph_t) = cph_counts(c);
    let (cpg_u, cpg_t) = (c.unconv_for(Context::CpG), c.total_for(Context::CpG));
    let readout = chemistry.readout();
    let calls_unconverted = readout.calls_unconverted();
    SummaryRow {
        sample: sample.to_string(),
        methylsieve_version: METHYLSIEVE_BUILD,
        chemistry: chemistry.name(),
        scope: scope.name.clone(),
        r1_mask_5p: masks.r1_5p,
        r2_mask_5p: masks.r2_5p,
        se_mask_5p: masks.se_5p,
        se_mask_3p: masks.se_3p,
        n_templates: scope.n_templates,
        n_mapped: scope.n_mapped,
        n_evaluated: scope.n_evaluated,
        n_unconverted: calls_unconverted.then_some(scope.n_unconverted),
        frac_unconverted: (calls_unconverted && scope.n_evaluated > 0)
            .then(|| scope.n_unconverted as f64 / scope.n_evaluated as f64),
        chimeric_to_control_templates: chimeric,
        cpa_obs: c.total_for(Context::CpA),
        cpa_conv_rate: conv_rate(c.unconv_for(Context::CpA), c.total_for(Context::CpA)),
        cpc_obs: c.total_for(Context::CpC),
        cpc_conv_rate: conv_rate(c.unconv_for(Context::CpC), c.total_for(Context::CpC)),
        cpt_obs: c.total_for(Context::CpT),
        cpt_conv_rate: conv_rate(c.unconv_for(Context::CpT), c.total_for(Context::CpT)),
        cph_obs: cph_t,
        cph_conv_rate: conv_rate(cph_u, cph_t),
        cpg_obs: cpg_t,
        cpg_conv_rate: conv_rate(cpg_u, cpg_t),
        cpg_meth_rate: (cpg_t > 0).then(|| readout.methylated(cpg_u, cpg_t) as f64 / cpg_t as f64),
    }
}

/// Build the summary rows: the genome scope first, then each control contig.
/// The applied mask plan rides along as run-level columns on every row.
fn summary_rows(
    stats: &Stats,
    chemistry: Chemistry,
    mask_plan: Option<&MaskPlan>,
    sample: &str,
) -> Vec<SummaryRow> {
    let masks = MaskCols::from_plan(mask_plan);
    let mut rows = vec![summary_row(
        sample,
        chemistry,
        &stats.genome,
        Some(stats.chimeric_to_control_templates),
        masks,
    )];
    rows.extend(stats.controls.iter().map(|c| summary_row(sample, chemistry, c, None, masks)));
    rows
}

// ── conversion-matrix.tsv ─────────────────────────────────────────────────────

/// One conversion-matrix cell: a `(checked, converted)` count of templates in one
/// scope, with the verdict the decision engine would assign that cell. The
/// verdict is blank for chemistries that make no unconverted call.
#[derive(Serialize)]
struct MatrixRow {
    sample: String,
    scope: String,
    context: String,
    checked_sites: u64,
    converted_sites: u64,
    #[serde(serialize_with = "ser_opt_f64_6dp")]
    conversion_rate: Option<f64>,
    n_templates: u64,
    decision: Option<&'static str>,
    decided_by: Option<&'static str>,
}

/// The scopes with a conversion matrix, in report order: the genome when the
/// chemistry makes the unconverted call, then every control. Each page carries
/// the boundary only when that call is made.
fn matrix_pages<'a>(
    stats: &'a Stats,
    calls_unconverted: bool,
    contexts: &'a str,
    control_contexts: &'a str,
) -> Vec<MatrixPage<'a>> {
    let genome = calls_unconverted.then(|| MatrixPage {
        scope: &stats.genome.name,
        is_control: false,
        contexts,
        cells: &stats.genome.conversion_matrix,
        boundary: true,
    });
    let controls = stats.controls.iter().map(|c| MatrixPage {
        scope: &c.name,
        is_control: true,
        contexts: control_contexts,
        cells: &c.conversion_matrix,
        boundary: calls_unconverted,
    });
    genome.into_iter().chain(controls).collect()
}

fn matrix_rows<F>(pages: &[MatrixPage<'_>], sample: &str, classify: &F) -> Vec<MatrixRow>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    let mut rows = Vec::new();
    for page in pages {
        for (&(checked, unconv), &n_templates) in page.cells {
            let verdict = page.boundary.then(|| classify(unconv, checked));
            let converted = checked - unconv;
            rows.push(MatrixRow {
                sample: sample.to_string(),
                scope: page.scope.to_string(),
                context: page.contexts.to_string(),
                checked_sites: checked,
                converted_sites: converted,
                conversion_rate: (checked > 0).then(|| converted as f64 / checked as f64),
                n_templates,
                decision: verdict.map(|(u, _)| if u { "unconverted" } else { "converted" }),
                decided_by: verdict.map(|(_, by)| by.as_str()),
            });
        }
    }
    rows
}

// ── mbias.tsv ─────────────────────────────────────────────────────────────────

/// One per-read-cycle methylation point, read through the chemistry. `ci_lo`/
/// `ci_hi` are a 95% Agresti–Coull interval on `frac_methylation`.
#[derive(Serialize)]
struct MbiasRow {
    sample: String,
    read: &'static str,
    end: &'static str,
    context: &'static str,
    cycle: usize,
    n_methylated: u64,
    n_total: u64,
    #[serde(serialize_with = "ser_f64_6dp")]
    frac_methylation: f64,
    #[serde(serialize_with = "ser_f64_6dp")]
    ci_lo: f64,
    #[serde(serialize_with = "ser_f64_6dp")]
    ci_hi: f64,
}

/// Whether `(role, end)` is a curve we report: paired/orphan reads only have a
/// meaningful 5' M-bias here (their 3' end is the mate's domain); single-end
/// reads track both ends.
fn reported(role: ReadRole, end: ReadEnd) -> bool {
    end == ReadEnd::FivePrime || role == ReadRole::Se
}

fn mbias_rows(mbias: &MbiasAccumulator, sample: &str) -> Vec<MbiasRow> {
    let readout = mbias.readout();
    let mut rows = Vec::new();
    for &role in &ReadRole::ALL {
        for &end in &ReadEnd::ALL {
            if !reported(role, end) || !mbias.has_data(role, end) {
                continue;
            }
            for &ctx in &Context::ALL {
                for (cycle, cc) in mbias.cycles(role, end, ctx).iter().enumerate() {
                    let Some(frac) = cc.frac_methylated(readout) else { continue };
                    let methylated = cc.methylated(readout);
                    let (lo, hi) = agresti_coull(methylated, cc.total());
                    rows.push(MbiasRow {
                        sample: sample.to_string(),
                        read: role.label(),
                        end: end.label(),
                        context: ctx.label(),
                        cycle: cycle + 1,
                        n_methylated: methylated,
                        n_total: cc.total(),
                        frac_methylation: frac,
                        ci_lo: lo,
                        ci_hi: hi,
                    });
                }
            }
        }
    }
    rows
}

/// 95% Agresti–Coull confidence interval for a binomial fraction `x/n`.
fn agresti_coull(x: u64, n: u64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 0.0);
    }
    const Z: f64 = 1.96;
    let n_t = n as f64 + Z * Z;
    let p_t = (x as f64 + Z * Z / 2.0) / n_t;
    let margin = Z * (p_t * (1.0 - p_t) / n_t).sqrt();
    ((p_t - margin).max(0.0), (p_t + margin).min(1.0))
}

// ── sample resolution ─────────────────────────────────────────────────────────

/// Resolve the `sample` value, in precedence order: explicit `--sample` →
/// comma-joined unique `@RG SM:` → the input file's stem (only when reading a
/// regular file, not stdin/a pipe) → `"unknown"`.
fn resolve_sample(
    header: &Header,
    sample_override: Option<&str>,
    input_path: Option<&Path>,
) -> String {
    if let Some(s) = sample_override {
        return s.to_string();
    }
    let mut samples: BTreeSet<String> = BTreeSet::new();
    for (_id, map) in header.read_groups() {
        if let Some(sm) = map.other_fields().get(&rg_tag::SAMPLE) {
            let s = sm.to_string();
            if !s.is_empty() {
                samples.insert(s);
            }
        }
    }
    if !samples.is_empty() {
        return samples.into_iter().collect::<Vec<_>>().join(",");
    }
    // Fall back to the input file stem, but only for a real regular file —
    // stdin/pipes/devices have no meaningful name.
    if let Some(stem) = input_path
        .filter(|p| std::fs::metadata(p).is_ok_and(|m| m.is_file()))
        .and_then(Path::file_stem)
        .and_then(|s| s.to_str())
        .filter(|s| !s.is_empty())
    {
        return stem.to_string();
    }
    "unknown".to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::Readout;
    use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
    use crate::sieve::{DecidedBy, PerContextCounters, Stats};

    fn genome_with(counters: PerContextCounters, n_eval: u64, n_unconv: u64) -> Stats {
        let mut s = Stats::new(&[]);
        s.genome.n_templates = n_eval;
        s.genome.n_mapped = n_eval;
        s.genome.n_evaluated = n_eval;
        s.genome.n_unconverted = n_unconv;
        s.genome.counters = counters;
        s
    }

    /// Serialize rows through the real TSV writer and return the file text, so
    /// tests assert on actual headers, ordering, and number formatting.
    fn render<T: Serialize>(rows: &[T]) -> String {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        DelimFile::default().write_tsv(tmp.path(), rows.iter()).unwrap();
        std::fs::read_to_string(tmp.path()).unwrap()
    }

    /// Column index of a header name in a rendered TSV's header line.
    fn col(header: &str, name: &str) -> usize {
        header.split('\t').position(|c| c == name).unwrap_or_else(|| panic!("no column {name}"))
    }

    #[test]
    fn summary_is_one_decision_row_per_scope() {
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpA, 1, 1000); // conv rate 0.999
        let stats = genome_with(counters, 100, 3);
        let text = render(&summary_rows(&stats, Chemistry::EmSeq, None, "S1"));
        let mut lines = text.lines();
        let hdr = lines.next().unwrap();
        // Summary is one decision-basis row per scope; there is no `read` column.
        assert!(!hdr.split('\t').any(|c| c == "read"), "no per-read column");
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "scope")], "genome");
        assert_eq!(row[col(hdr, "n_templates")], "100");
        assert_eq!(row[col(hdr, "n_unconverted")], "3");
        assert_eq!(row[col(hdr, "frac_unconverted")], "0.030000");
        assert_eq!(row[col(hdr, "CpA_obs")], "1000");
        assert_eq!(row[col(hdr, "CpA_conv_rate")], "0.999000");
        // No mask plan → all four mask columns blank.
        assert_eq!(row[col(hdr, "r1_mask_5p")], "");
        assert_eq!(row[col(hdr, "se_mask_3p")], "");
        // Header + the single genome row only (no controls configured).
        assert_eq!(text.lines().count(), 2);
    }

    #[test]
    fn mask_plan_fills_mask_columns_on_every_scope_row() {
        let mut stats = Stats::new(&["chrCtrl".to_string()]);
        stats.genome.n_templates = 10;
        let plan = MaskPlan::explicit(2, 22, 2);
        let text = render(&summary_rows(&stats, Chemistry::EmSeq, Some(&plan), "S1"));
        let hdr = text.lines().next().unwrap().to_string();
        for row in text.lines().skip(1) {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols[col(&hdr, "r1_mask_5p")], "2");
            assert_eq!(cols[col(&hdr, "r2_mask_5p")], "22");
            // SE lengths default to 0 in an explicit (paired) plan.
            assert_eq!(cols[col(&hdr, "se_mask_5p")], "0");
            assert_eq!(cols[col(&hdr, "se_mask_3p")], "0");
        }
        // Header + genome + one control row.
        assert_eq!(text.lines().count(), 3);
    }

    #[test]
    fn mbias_tsv_has_per_cycle_rows() {
        let mut mbias = MbiasAccumulator::new(Readout::UnmethylatedConverts);
        for i in 0..100u64 {
            mbias.record(ReadRole::R2, ReadEnd::FivePrime, Context::CpG, 0, i < 30);
        }
        let text = render(&mbias_rows(&mbias, "S1"));
        let hdr = text.lines().next().unwrap();
        let row = text.lines().find(|l| l.contains("\tR2\t5p\tCpG\t1\t")).unwrap();
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[col(hdr, "n_methylated")], "30");
        assert_eq!(cols[col(hdr, "n_total")], "100");
        assert_eq!(cols[col(hdr, "frac_methylation")], "0.300000");
    }

    /// The decision used by the matrix tests: proportion at/above 40 sites, count below.
    fn test_classify(unconv: u64, checked: u64) -> (bool, DecidedBy) {
        if checked == 0 {
            (false, DecidedBy::TooFewSites)
        } else if checked >= 40 {
            (unconv as f64 / checked as f64 > 0.05, DecidedBy::Proportion)
        } else {
            (unconv >= 3, DecidedBy::Count)
        }
    }

    #[test]
    fn matrix_renders_sorted_cells() {
        let mut stats = Stats::new(&[]);
        stats.genome.conversion_matrix.insert((50, 5), 2);
        stats.genome.conversion_matrix.insert((0, 0), 5);
        stats.genome.conversion_matrix.insert((10, 3), 7);
        let pages = matrix_pages(&stats, true, "CpH", "CpH");
        let text = render(&matrix_rows(&pages, "S1", &test_classify));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "sample\tscope\tcontext\tchecked_sites\tconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
        );
        // converted_sites = checked - unconverted; conversion_rate = converted/checked.
        assert_eq!(lines[1], "S1\tgenome\tCpH\t0\t0\t\t5\tconverted\ttoo_few_sites");
        assert_eq!(lines[2], "S1\tgenome\tCpH\t10\t7\t0.700000\t7\tunconverted\tcount");
        assert_eq!(lines[3], "S1\tgenome\tCpH\t50\t45\t0.900000\t2\tunconverted\tproportion");
    }

    #[test]
    fn emseq_matrix_gives_controls_their_own_rows_and_verdicts() {
        let mut stats = Stats::new(&["lambda".to_string()]);
        stats.genome.conversion_matrix.insert((10, 0), 4);
        stats.controls[0].conversion_matrix.insert((10, 3), 1);
        let pages = matrix_pages(&stats, true, "CpH", "CpH");
        let text = render(&matrix_rows(&pages, "S1", &test_classify));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[1], "S1\tgenome\tCpH\t10\t10\t1.000000\t4\tconverted\tcount");
        assert_eq!(lines[2], "S1\tlambda\tCpH\t10\t7\t0.700000\t1\tunconverted\tcount");
    }

    #[test]
    fn taps_matrix_has_only_control_rows_with_blank_verdicts() {
        let mut stats = Stats::new(&["lambda".to_string()]);
        stats.genome.conversion_matrix.insert((10, 0), 4);
        stats.controls[0].conversion_matrix.insert((20, 18), 1);
        let pages = matrix_pages(&stats, false, "CpH", "CpG");
        let text = render(&matrix_rows(&pages, "S1", &test_classify));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1], "S1\tlambda\tCpG\t20\t2\t0.100000\t1\t\t");
    }

    #[test]
    fn taps_without_controls_has_no_matrix() {
        let stats = Stats::new(&[]);
        assert!(matrix_pages(&stats, false, "CpH", "CpG").is_empty());
    }

    #[test]
    fn emseq_summary_counts_unconverted_control_templates() {
        let mut stats = Stats::new(&["lambda".to_string()]);
        stats.controls[0].n_evaluated = 10;
        stats.controls[0].n_unconverted = 2;
        let text = render(&summary_rows(&stats, Chemistry::EmSeq, None, "S1"));
        let hdr = text.lines().next().unwrap();
        let row: Vec<&str> = text.lines().nth(2).unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "scope")], "lambda");
        assert_eq!(row[col(hdr, "n_unconverted")], "2");
        assert_eq!(row[col(hdr, "frac_unconverted")], "0.200000");
    }

    #[test]
    fn summary_records_the_chemistry() {
        let stats = genome_with(PerContextCounters::default(), 1, 0);
        let text = render(&summary_rows(&stats, Chemistry::Illumina5Base, None, "S1"));
        let hdr = text.lines().next().unwrap();
        let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "chemistry")], "illumina-5base");
    }

    #[test]
    fn taps_summary_leaves_decision_columns_blank() {
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpG, 40, 100);
        let stats = genome_with(counters, 10, 0);
        let text = render(&summary_rows(&stats, Chemistry::Taps, None, "S1"));
        let hdr = text.lines().next().unwrap();
        let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "n_unconverted")], "");
        assert_eq!(row[col(hdr, "frac_unconverted")], "");
        assert_eq!(row[col(hdr, "n_evaluated")], "10");
    }

    #[test]
    fn taps_summary_cpg_meth_rate_equals_cpg_conv_rate() {
        // 40 of 100 CpGs retained: conversion 0.6, which TAPS reads as methylation.
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpG, 40, 100);
        let stats = genome_with(counters, 10, 0);
        let text = render(&summary_rows(&stats, Chemistry::Taps, None, "S1"));
        let hdr = text.lines().next().unwrap();
        let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "CpG_conv_rate")], "0.600000");
        assert_eq!(row[col(hdr, "CpG_meth_rate")], "0.600000");
    }

    #[test]
    fn emseq_summary_cpg_meth_rate_is_the_retained_fraction() {
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpG, 40, 100);
        let stats = genome_with(counters, 10, 0);
        let text = render(&summary_rows(&stats, Chemistry::EmSeq, None, "S1"));
        let hdr = text.lines().next().unwrap();
        let row: Vec<&str> = text.lines().nth(1).unwrap().split('\t').collect();
        assert_eq!(row[col(hdr, "CpG_meth_rate")], "0.400000");
    }

    #[test]
    fn taps_mbias_tsv_counts_converted_cytosines_as_methylated() {
        let mut mbias = MbiasAccumulator::new(Readout::MethylatedConverts);
        for i in 0..100u64 {
            mbias.record(ReadRole::R2, ReadEnd::FivePrime, Context::CpG, 0, i < 30);
        }
        let text = render(&mbias_rows(&mbias, "S1"));
        let hdr = text.lines().next().unwrap();
        let row = text.lines().find(|l| l.contains("\tR2\t5p\tCpG\t1\t")).unwrap();
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[col(hdr, "n_methylated")], "70");
        assert_eq!(cols[col(hdr, "frac_methylation")], "0.700000");
    }
}
