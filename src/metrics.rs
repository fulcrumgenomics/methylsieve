//! Metric TSVs written under `--metrics-prefix PREFIX`:
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
//!
//! All rates are fractions in `[0, 1]` (never percentages).

use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context as _, Result};
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

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
/// accumulated during the run; `classify(unconv, monitored)` replays the
/// per-cell decision verdict for the conversion matrix.
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
    write_file(&with_suffix(prefix, "summary.tsv"), |w| {
        write_summary(w, stats, mask_plan, &sample)
    })?;
    write_file(&with_suffix(prefix, "conversion-matrix.tsv"), |w| {
        write_matrix(w, stats, &sample, &classify)
    })?;
    write_file(&with_suffix(prefix, "mbias.tsv"), |w| write_mbias(w, mbias, &sample))?;
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

/// Create `path` and run `render` against it, with path context on error.
fn write_file<F>(path: &Path, render: F) -> Result<()>
where
    F: FnOnce(&mut dyn Write) -> std::io::Result<()>,
{
    let mut f =
        std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    render(&mut f).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

// ── summary.tsv ───────────────────────────────────────────────────────────────

/// One summary row — one per scope, all decision-basis.
struct Row<'a> {
    sample: &'a str,
    scope: &'a str,
    n_templates: u64,
    n_mapped: u64,
    n_evaluated: u64,
    n_unconverted: u64,
    /// Chimeric-to-control count — only the genome scope (control diagnostic);
    /// `None` for control scopes.
    chimeric_to_control: Option<u64>,
    counters: PerContextCounters,
    /// Applied per-read mask lengths from the frozen plan (a run-level property
    /// repeated on every scope row), blank when masking was not run.
    masks: MaskCols,
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

type ColumnFn = fn(&Row) -> String;

/// Optional `u64` → string, blank when `None`.
fn opt(v: Option<u64>) -> String {
    v.map(|n| n.to_string()).unwrap_or_default()
}

/// Optional `usize` → string, blank when `None`.
fn opt_usize(v: Option<usize>) -> String {
    v.map(|n| n.to_string()).unwrap_or_default()
}

/// Conversion rate `1 - unconv/total` from explicit counts, blank when no sites.
fn conv_rate(unconv: u64, total: u64) -> String {
    if total == 0 { String::new() } else { format!("{:.6}", 1.0 - unconv as f64 / total as f64) }
}

/// Combined CpH `(unconv, total)` across CpA/CpC/CpT.
fn cph_counts(c: &PerContextCounters) -> (u64, u64) {
    [Context::CpA, Context::CpC, Context::CpT]
        .iter()
        .fold((0, 0), |(u, t), &ctx| (u + c.unconv_for(ctx), t + c.total_for(ctx)))
}

const COLUMNS: &[(&str, ColumnFn)] = &[
    ("sample", |r| r.sample.to_string()),
    ("methylsieve_version", |_| METHYLSIEVE_BUILD.to_string()),
    ("scope", |r| r.scope.to_string()),
    ("r1_mask_5p", |r| opt_usize(r.masks.r1_5p)),
    ("r2_mask_5p", |r| opt_usize(r.masks.r2_5p)),
    ("se_mask_5p", |r| opt_usize(r.masks.se_5p)),
    ("se_mask_3p", |r| opt_usize(r.masks.se_3p)),
    ("n_templates", |r| r.n_templates.to_string()),
    ("n_mapped", |r| r.n_mapped.to_string()),
    ("n_evaluated", |r| r.n_evaluated.to_string()),
    ("n_unconverted", |r| r.n_unconverted.to_string()),
    ("frac_unconverted", |r| {
        if r.n_evaluated > 0 {
            format!("{:.6}", r.n_unconverted as f64 / r.n_evaluated as f64)
        } else {
            String::new()
        }
    }),
    ("chimeric_to_control_templates", |r| opt(r.chimeric_to_control)),
    ("CpA_obs", |r| r.counters.total_for(Context::CpA).to_string()),
    ("CpA_conv_rate", |r| {
        conv_rate(r.counters.unconv_for(Context::CpA), r.counters.total_for(Context::CpA))
    }),
    ("CpC_obs", |r| r.counters.total_for(Context::CpC).to_string()),
    ("CpC_conv_rate", |r| {
        conv_rate(r.counters.unconv_for(Context::CpC), r.counters.total_for(Context::CpC))
    }),
    ("CpT_obs", |r| r.counters.total_for(Context::CpT).to_string()),
    ("CpT_conv_rate", |r| {
        conv_rate(r.counters.unconv_for(Context::CpT), r.counters.total_for(Context::CpT))
    }),
    ("CpH_obs", |r| cph_counts(&r.counters).1.to_string()),
    ("CpH_conv_rate", |r| {
        let (u, t) = cph_counts(&r.counters);
        conv_rate(u, t)
    }),
    ("CpG_obs", |r| r.counters.total_for(Context::CpG).to_string()),
    ("CpG_conv_rate", |r| {
        conv_rate(r.counters.unconv_for(Context::CpG), r.counters.total_for(Context::CpG))
    }),
    // CpG methylation rate = unconverted/total CpG = 1 - CpG conv rate, stated
    // directly because it is the headline biological readout.
    ("CpG_meth_rate", |r| {
        let total = r.counters.total_for(Context::CpG);
        if total == 0 {
            String::new()
        } else {
            format!("{:.6}", r.counters.unconv_for(Context::CpG) as f64 / total as f64)
        }
    }),
];

/// One scope row (decision basis). `chimeric` is the genome-only control
/// diagnostic (`None` for control scopes); `masks` is the run-level mask plan,
/// repeated on every row.
fn scope_row<'a>(
    sample: &'a str,
    scope: &'a ScopeStats,
    chimeric: Option<u64>,
    masks: MaskCols,
) -> Row<'a> {
    Row {
        sample,
        scope: &scope.name,
        n_templates: scope.n_templates,
        n_mapped: scope.n_mapped,
        n_evaluated: scope.n_evaluated,
        n_unconverted: scope.n_unconverted,
        chimeric_to_control: chimeric,
        counters: scope.counters,
        masks,
    }
}

/// Render the summary: a header row, then one decision-basis row per scope
/// (genome first, then each control contig). The applied mask plan rides along
/// as run-level columns on every row.
fn write_summary(
    w: &mut dyn Write,
    stats: &Stats,
    mask_plan: Option<&MaskPlan>,
    sample: &str,
) -> std::io::Result<()> {
    let header: String = COLUMNS.iter().map(|(n, _)| *n).collect::<Vec<_>>().join("\t");
    writeln!(w, "{header}")?;

    let masks = MaskCols::from_plan(mask_plan);
    let render = |w: &mut dyn Write, row: &Row| -> std::io::Result<()> {
        let line: String = COLUMNS.iter().map(|(_, f)| f(row)).collect::<Vec<_>>().join("\t");
        writeln!(w, "{line}")
    };

    render(w, &scope_row(sample, &stats.genome, Some(stats.chimeric_to_control_templates), masks))?;
    for control in &stats.controls {
        render(w, &scope_row(sample, control, None, masks))?;
    }
    Ok(())
}

// ── conversion-matrix.tsv ─────────────────────────────────────────────────────

fn write_matrix<F>(
    w: &mut dyn Write,
    stats: &Stats,
    sample: &str,
    classify: &F,
) -> std::io::Result<()>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    writeln!(
        w,
        "sample\tchecked_sites\tconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
    )?;
    for (&(checked, unconv), &n_templates) in &stats.conversion_matrix {
        let (unconverted, by) = classify(unconv, checked);
        let converted = checked - unconv;
        let conv_rate = if checked > 0 {
            format!("{:.6}", converted as f64 / checked as f64)
        } else {
            String::new()
        };
        let decision = if unconverted { "unconverted" } else { "converted" };
        writeln!(
            w,
            "{sample}\t{checked}\t{converted}\t{conv_rate}\t{n_templates}\t{decision}\t{}",
            by.as_str()
        )?;
    }
    Ok(())
}

// ── mbias.tsv ─────────────────────────────────────────────────────────────────

/// Whether `(role, end)` is a curve we report: paired/orphan reads only have a
/// meaningful 5' M-bias here (their 3' end is the mate's domain); single-end
/// reads track both ends.
fn reported(role: ReadRole, end: ReadEnd) -> bool {
    end == ReadEnd::FivePrime || role == ReadRole::Se
}

/// Per-read-cycle methylation curve: one row per `(read, end, context, cycle)`
/// with coverage. `ci_lo`/`ci_hi` are a 95% Agresti–Coull interval.
fn write_mbias(w: &mut dyn Write, mbias: &MbiasAccumulator, sample: &str) -> std::io::Result<()> {
    writeln!(
        w,
        "sample\tread\tend\tcontext\tcycle\tn_methylated\tn_total\tfrac_methylation\tci_lo\tci_hi"
    )?;
    for &role in &ReadRole::ALL {
        for &end in &ReadEnd::ALL {
            if !reported(role, end) || !mbias.has_data(role, end) {
                continue;
            }
            for &ctx in &Context::ALL {
                for (cycle, cc) in mbias.cycles(role, end, ctx).iter().enumerate() {
                    let Some(frac) = cc.frac() else { continue };
                    let (lo, hi) = agresti_coull(cc.meth(), cc.total());
                    writeln!(
                        w,
                        "{sample}\t{}\t{}\t{}\t{}\t{}\t{}\t{frac:.6}\t{lo:.6}\t{hi:.6}",
                        role.label(),
                        end.label(),
                        ctx.label(),
                        cycle + 1,
                        cc.meth(),
                        cc.total(),
                    )?;
                }
            }
        }
    }
    Ok(())
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

    fn render_summary(stats: &Stats, mask_plan: Option<&MaskPlan>) -> String {
        let mut buf = Vec::new();
        write_summary(&mut buf, stats, mask_plan, "S1").unwrap();
        String::from_utf8(buf).unwrap()
    }

    /// Column index of a header name in the summary schema.
    fn col(name: &str) -> usize {
        COLUMNS.iter().position(|(n, _)| *n == name).unwrap_or_else(|| panic!("no column {name}"))
    }

    #[test]
    fn summary_is_one_decision_row_per_scope() {
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpA, 1, 1000); // conv rate 0.999
        let stats = genome_with(counters, 100, 3);
        let text = render_summary(&stats, None);
        let mut lines = text.lines();
        let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
        // The per-read fold is gone — no `read` column, one row per scope.
        assert!(!hdr.contains(&"read"), "read column removed");
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(hdr.len(), row.len());
        assert_eq!(row[col("scope")], "genome");
        assert_eq!(row[col("n_templates")], "100");
        assert_eq!(row[col("n_unconverted")], "3");
        assert_eq!(row[col("frac_unconverted")], "0.030000");
        assert_eq!(row[col("CpA_obs")], "1000");
        assert_eq!(row[col("CpA_conv_rate")], "0.999000");
        // No mask plan → all four mask columns blank.
        assert_eq!(row[col("r1_mask_5p")], "");
        assert_eq!(row[col("se_mask_3p")], "");
        // Header + the single genome row only (no controls configured).
        assert_eq!(text.lines().count(), 2);
    }

    #[test]
    fn mask_plan_fills_mask_columns_on_every_scope_row() {
        let mut stats = Stats::new(&["chrCtrl".to_string()]);
        stats.genome.n_templates = 10;
        let plan = MaskPlan::explicit(2, 22, 2);
        let text = render_summary(&stats, Some(&plan));
        for row in text.lines().skip(1) {
            let cols: Vec<&str> = row.split('\t').collect();
            assert_eq!(cols[col("r1_mask_5p")], "2");
            assert_eq!(cols[col("r2_mask_5p")], "22");
            // SE lengths default to 0 in an explicit (paired) plan.
            assert_eq!(cols[col("se_mask_5p")], "0");
            assert_eq!(cols[col("se_mask_3p")], "0");
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
        let mut buf = Vec::new();
        write_mbias(&mut buf, &mbias, "S1").unwrap();
        let text = String::from_utf8(buf).unwrap();
        let row = text.lines().find(|l| l.contains("\tR2\t5p\tCpG\t1\t")).unwrap();
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[5], "30");
        assert_eq!(cols[6], "100");
        assert_eq!(cols[7], "0.300000");
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
        let mut buf = Vec::new();
        write_matrix(&mut buf, &stats, "S1", &classify).unwrap();
        let lines: Vec<String> =
            String::from_utf8(buf).unwrap().lines().map(str::to_string).collect();
        // converted_sites = checked - unconverted; conversion_rate = converted/checked.
        assert_eq!(lines[1], "S1\t0\t0\t\t5\tconverted\ttoo_few_sites");
        assert_eq!(lines[2], "S1\t10\t7\t0.700000\t7\tunconverted\tcount");
        assert_eq!(lines[3], "S1\t50\t45\t0.900000\t2\tunconverted\tproportion");
    }
}
