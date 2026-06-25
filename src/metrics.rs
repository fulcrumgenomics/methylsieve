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
//! - `PREFIX.conversion-matrix.tsv` — per-`(checked, converted)` decision cell.
//! - `PREFIX.mbias.tsv` — per-read-cycle methylation by `(read, end, context)`.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use anyhow::{Result, anyhow};
use fgoxide::io::DelimFile;
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;
use serde::Serialize;

use crate::METHYLSIEVE_BUILD;
use crate::mask::MaskPlan;
use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
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
/// decision verdict for the conversion matrix.
#[allow(clippy::too_many_arguments)] // a single wiring call site; a params struct would not aid clarity
pub(crate) fn write_all<F>(
    prefix: &Path,
    stats: &Stats,
    mbias: &MbiasAccumulator,
    mask_plan: Option<&MaskPlan>,
    header: &Header,
    sample_override: Option<&str>,
    input_path: Option<&Path>,
    contexts: &str,
    classify: F,
) -> Result<()>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    let sample = resolve_sample(header, sample_override, input_path);
    let df = DelimFile::default();
    write_tsv(&df, &with_suffix(prefix, "summary.tsv"), summary_rows(stats, mask_plan, &sample))?;
    write_tsv(
        &df,
        &with_suffix(prefix, "conversion-matrix.tsv"),
        matrix_rows(stats, &sample, &classify),
    )?;
    write_tsv(&df, &with_suffix(prefix, "mbias.tsv"), mbias_rows(mbias, &sample))?;
    // PDF plots of the same data (M-bias curves + conversion-matrix hexbin).
    crate::plots::write_mbias_pdf(&with_suffix(prefix, "mbias.pdf"), mbias, mask_plan, &sample)?;
    crate::plots::write_matrix_pdf(
        &with_suffix(prefix, "conversion-matrix.pdf"),
        &stats.conversion_matrix,
        &classify,
        &sample,
        contexts,
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
#[derive(Serialize)]
struct SummaryRow {
    sample: String,
    methylsieve_version: &'static str,
    scope: String,
    r1_mask_5p: Option<usize>,
    r2_mask_5p: Option<usize>,
    se_mask_5p: Option<usize>,
    se_mask_3p: Option<usize>,
    n_templates: u64,
    n_mapped: u64,
    n_evaluated: u64,
    n_unconverted: u64,
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
    // CpG methylation rate = unconverted/total CpG = 1 - CpG conv rate, stated
    // directly because it is the headline biological readout.
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
    scope: &ScopeStats,
    chimeric: Option<u64>,
    masks: MaskCols,
) -> SummaryRow {
    let c = &scope.counters;
    let (cph_u, cph_t) = cph_counts(c);
    let (cpg_u, cpg_t) = (c.unconv_for(Context::CpG), c.total_for(Context::CpG));
    SummaryRow {
        sample: sample.to_string(),
        methylsieve_version: METHYLSIEVE_BUILD,
        scope: scope.name.clone(),
        r1_mask_5p: masks.r1_5p,
        r2_mask_5p: masks.r2_5p,
        se_mask_5p: masks.se_5p,
        se_mask_3p: masks.se_3p,
        n_templates: scope.n_templates,
        n_mapped: scope.n_mapped,
        n_evaluated: scope.n_evaluated,
        n_unconverted: scope.n_unconverted,
        frac_unconverted: (scope.n_evaluated > 0)
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
        cpg_meth_rate: (cpg_t > 0).then(|| cpg_u as f64 / cpg_t as f64),
    }
}

/// Build the summary rows: the genome scope first, then each control contig.
/// The applied mask plan rides along as run-level columns on every row.
fn summary_rows(stats: &Stats, mask_plan: Option<&MaskPlan>, sample: &str) -> Vec<SummaryRow> {
    let masks = MaskCols::from_plan(mask_plan);
    let mut rows =
        vec![summary_row(sample, &stats.genome, Some(stats.chimeric_to_control_templates), masks)];
    rows.extend(stats.controls.iter().map(|c| summary_row(sample, c, None, masks)));
    rows
}

// ── conversion-matrix.tsv ─────────────────────────────────────────────────────

/// One conversion-matrix cell: a `(checked, converted)` count of templates with
/// the verdict the decision engine would assign that cell.
#[derive(Serialize)]
struct MatrixRow {
    sample: String,
    checked_sites: u64,
    converted_sites: u64,
    #[serde(serialize_with = "ser_opt_f64_6dp")]
    conversion_rate: Option<f64>,
    n_templates: u64,
    decision: &'static str,
    decided_by: &'static str,
}

fn matrix_rows<F>(stats: &Stats, sample: &str, classify: &F) -> Vec<MatrixRow>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    stats
        .conversion_matrix
        .iter()
        .map(|(&(checked, unconv), &n_templates)| {
            let (unconverted, by) = classify(unconv, checked);
            let converted = checked - unconv;
            MatrixRow {
                sample: sample.to_string(),
                checked_sites: checked,
                converted_sites: converted,
                conversion_rate: (checked > 0).then(|| converted as f64 / checked as f64),
                n_templates,
                decision: if unconverted { "unconverted" } else { "converted" },
                decided_by: by.as_str(),
            }
        })
        .collect()
}

// ── mbias.tsv ─────────────────────────────────────────────────────────────────

/// One per-read-cycle methylation point. `ci_lo`/`ci_hi` are a 95%
/// Agresti–Coull interval on `frac_methylation`.
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
    let mut rows = Vec::new();
    for &role in &ReadRole::ALL {
        for &end in &ReadEnd::ALL {
            if !reported(role, end) || !mbias.has_data(role, end) {
                continue;
            }
            for &ctx in &Context::ALL {
                for (cycle, cc) in mbias.cycles(role, end, ctx).iter().enumerate() {
                    let Some(frac) = cc.frac() else { continue };
                    let (lo, hi) = agresti_coull(cc.meth(), cc.total());
                    rows.push(MbiasRow {
                        sample: sample.to_string(),
                        read: role.label(),
                        end: end.label(),
                        context: ctx.label(),
                        cycle: cycle + 1,
                        n_methylated: cc.meth(),
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
        let text = render(&summary_rows(&stats, None, "S1"));
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
        let text = render(&summary_rows(&stats, Some(&plan), "S1"));
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
        let mut mbias = MbiasAccumulator::new();
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

    #[test]
    fn matrix_renders_sorted_cells() {
        let mut stats = Stats::new(&[]);
        stats.conversion_matrix.insert((50, 5), 2);
        stats.conversion_matrix.insert((0, 0), 5);
        stats.conversion_matrix.insert((10, 3), 7);
        let classify = |unconv: u64, checked: u64| -> (bool, DecidedBy) {
            if checked == 0 {
                (false, DecidedBy::TooFewSites)
            } else if checked >= 40 {
                (unconv as f64 / checked as f64 > 0.05, DecidedBy::Proportion)
            } else {
                (unconv >= 3, DecidedBy::Count)
            }
        };
        let text = render(&matrix_rows(&stats, "S1", &classify));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "sample\tchecked_sites\tconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
        );
        // converted_sites = checked - unconverted; conversion_rate = converted/checked.
        assert_eq!(lines[1], "S1\t0\t0\t\t5\tconverted\ttoo_few_sites");
        assert_eq!(lines[2], "S1\t10\t7\t0.700000\t7\tunconverted\tcount");
        assert_eq!(lines[3], "S1\t50\t45\t0.900000\t2\tunconverted\tproportion");
    }
}
