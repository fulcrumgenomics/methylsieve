//! Metric TSVs written under `--metrics-prefix PREFIX`:
//!
//! - `PREFIX.summary.tsv` — per-context conversion summary, **folded** over a
//!   `read` dimension: one `all` row per scope (genome first, then each
//!   `--control-contig`) carrying the per-template decision counts and run
//!   diagnostics, followed by `R1`/`R2`/`SE` rows (genome only) with per-read
//!   conversion broken out. The `all` row is decision-basis (overlap-deduped,
//!   end-trimmed, includes supplementary evidence); the per-read rows are
//!   M-bias-basis (every primary call, base-quality-gated, no dedup), so they
//!   intentionally need not sum to `all`.
//! - `PREFIX.conversion_matrix.tsv` — per-`(checked, unconverted)` decision cell.
//! - `PREFIX.mbias.tsv` — per-read-cycle methylation by `(read, end, context)`.
//! - `PREFIX.mbias_bounds.tsv` — suggested 5'/3' mask lengths + plateau per read.
//!
//! All rates are fractions in `[0, 1]` (never percentages).

use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context as _, Result};
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

use crate::METHYLSIEVE_BUILD;
use crate::mbias::{DetectParams, MbiasAccumulator, ReadEnd, ReadRole, detect_mask_length};
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
pub(crate) fn write_all<F>(
    prefix: &Path,
    stats: &Stats,
    mbias: &MbiasAccumulator,
    detect: DetectParams,
    header: &Header,
    sample_override: Option<&str>,
    classify: F,
) -> Result<()>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    let sample = resolve_sample(header, sample_override);
    write_file(&with_suffix(prefix, "summary.tsv"), |w| write_summary(w, stats, mbias, &sample))?;
    write_file(&with_suffix(prefix, "conversion_matrix.tsv"), |w| {
        write_matrix(w, stats, &sample, &classify)
    })?;
    write_file(&with_suffix(prefix, "mbias.tsv"), |w| write_mbias(w, mbias, &sample))?;
    write_file(&with_suffix(prefix, "mbias_bounds.tsv"), |w| {
        write_mbias_bounds(w, mbias, detect, &sample)
    })?;
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

/// Whole-run diagnostic counters, rendered only on the genome `all` row.
#[derive(Debug, Clone, Copy)]
struct Diagnostics {
    chimeric_to_control_templates: u64,
    unmapped_templates: u64,
    zero_site_templates: u64,
    below_min_sites_templates: u64,
}

/// One summary row. Decision-level counts are `None` for per-read rows (those
/// quantities are per-template, not per-read).
struct Row<'a> {
    sample: &'a str,
    scope: &'a str,
    read: &'a str,
    n_templates: Option<u64>,
    n_evaluated: Option<u64>,
    n_unconverted: Option<u64>,
    n_removed: Option<u64>,
    counters: PerContextCounters,
    diagnostics: Option<Diagnostics>,
}

type ColumnFn = fn(&Row) -> String;

/// Optional `u64` → string, blank when `None`.
fn opt(v: Option<u64>) -> String {
    v.map(|n| n.to_string()).unwrap_or_default()
}

/// `1 - unconv/total` for `ctx`, blank when no sites.
fn conv_rate(c: &PerContextCounters, ctx: Context) -> String {
    let total = c.total_for(ctx);
    if total == 0 {
        String::new()
    } else {
        format!("{:.6}", 1.0 - c.unconv_for(ctx) as f64 / total as f64)
    }
}

const COLUMNS: &[(&str, ColumnFn)] = &[
    ("sample", |r| r.sample.to_string()),
    ("methylsieve_version", |_| METHYLSIEVE_BUILD.to_string()),
    ("scope", |r| r.scope.to_string()),
    ("read", |r| r.read.to_string()),
    ("n_templates", |r| opt(r.n_templates)),
    ("n_evaluated", |r| opt(r.n_evaluated)),
    ("n_unconverted", |r| opt(r.n_unconverted)),
    ("n_removed", |r| opt(r.n_removed)),
    ("frac_unconverted", |r| match (r.n_unconverted, r.n_evaluated) {
        (Some(u), Some(e)) if e > 0 => format!("{:.6}", u as f64 / e as f64),
        _ => String::new(),
    }),
    ("CA_unconv", |r| r.counters.unconv_for(Context::CpA).to_string()),
    ("CA_total", |r| r.counters.total_for(Context::CpA).to_string()),
    ("CC_unconv", |r| r.counters.unconv_for(Context::CpC).to_string()),
    ("CC_total", |r| r.counters.total_for(Context::CpC).to_string()),
    ("CT_unconv", |r| r.counters.unconv_for(Context::CpT).to_string()),
    ("CT_total", |r| r.counters.total_for(Context::CpT).to_string()),
    ("CG_unconv", |r| r.counters.unconv_for(Context::CpG).to_string()),
    ("CG_total", |r| r.counters.total_for(Context::CpG).to_string()),
    ("conv_rate_CpA", |r| conv_rate(&r.counters, Context::CpA)),
    ("conv_rate_CpC", |r| conv_rate(&r.counters, Context::CpC)),
    ("conv_rate_CpT", |r| conv_rate(&r.counters, Context::CpT)),
    ("conv_rate_CpG", |r| conv_rate(&r.counters, Context::CpG)),
    ("chimeric_to_control_templates", |r| {
        opt(r.diagnostics.map(|d| d.chimeric_to_control_templates))
    }),
    ("unmapped_templates", |r| opt(r.diagnostics.map(|d| d.unmapped_templates))),
    ("zero_site_templates", |r| opt(r.diagnostics.map(|d| d.zero_site_templates))),
    ("below_min_sites_templates", |r| opt(r.diagnostics.map(|d| d.below_min_sites_templates))),
];

/// The `all` row for a scope (decision basis).
fn all_row<'a>(
    sample: &'a str,
    scope: &'a ScopeStats,
    diagnostics: Option<Diagnostics>,
) -> Row<'a> {
    Row {
        sample,
        scope: &scope.name,
        read: "all",
        n_templates: Some(scope.n_templates),
        n_evaluated: Some(scope.n_evaluated),
        n_unconverted: Some(scope.n_unconverted),
        n_removed: Some(scope.n_removed),
        counters: scope.counters,
        diagnostics,
    }
}

/// Per-read context counters from the M-bias accumulator (5' end, summed over
/// all cycles): every primary call of that role, base-quality-gated. Summing the
/// 5' end alone counts each site exactly once — including single-end reads, where
/// every site is recorded under both ends, so the 5' sum is already the full
/// total (summing both ends would double-count).
fn role_counters(mbias: &MbiasAccumulator, role: ReadRole) -> PerContextCounters {
    let mut c = PerContextCounters::default();
    for &ctx in &Context::ALL {
        for cc in mbias.cycles(role, ReadEnd::FivePrime, ctx) {
            c.add_counts(ctx, cc.meth(), cc.total());
        }
    }
    c
}

/// Render the folded summary: a header row, then per scope an `all` row; the
/// genome scope additionally gets `R1`/`R2`/`SE` rows from the M-bias counts.
fn write_summary(
    w: &mut dyn Write,
    stats: &Stats,
    mbias: &MbiasAccumulator,
    sample: &str,
) -> std::io::Result<()> {
    let diagnostics = Diagnostics {
        chimeric_to_control_templates: stats.chimeric_to_control_templates,
        unmapped_templates: stats.unmapped_templates,
        zero_site_templates: stats.zero_site_templates,
        below_min_sites_templates: stats.below_min_sites_templates,
    };

    let header: String = COLUMNS.iter().map(|(n, _)| *n).collect::<Vec<_>>().join("\t");
    writeln!(w, "{header}")?;

    let render = |w: &mut dyn Write, row: &Row| -> std::io::Result<()> {
        let line: String = COLUMNS.iter().map(|(_, f)| f(row)).collect::<Vec<_>>().join("\t");
        writeln!(w, "{line}")
    };

    render(w, &all_row(sample, &stats.genome, Some(diagnostics)))?;
    // Per-read rows (genome only): M-bias basis. Emitted when the role has data.
    for &role in &ReadRole::ALL {
        if mbias.has_data(role, ReadEnd::FivePrime) {
            render(
                w,
                &Row {
                    sample,
                    scope: &stats.genome.name,
                    read: role.label(),
                    n_templates: None,
                    n_evaluated: None,
                    n_unconverted: None,
                    n_removed: None,
                    counters: role_counters(mbias, role),
                    diagnostics: None,
                },
            )?;
        }
    }
    for control in &stats.controls {
        render(w, &all_row(sample, control, None))?;
    }
    Ok(())
}

// ── conversion_matrix.tsv ─────────────────────────────────────────────────────

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
        "sample\tchecked_sites\tunconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
    )?;
    for (&(checked, unconv), &n_templates) in &stats.conversion_matrix {
        let (unconverted, by) = classify(unconv, checked);
        let conv_rate = if checked > 0 {
            format!("{:.6}", 1.0 - unconv as f64 / checked as f64)
        } else {
            String::new()
        };
        let decision = if unconverted { "unconverted" } else { "converted" };
        writeln!(
            w,
            "{sample}\t{checked}\t{unconv}\t{conv_rate}\t{n_templates}\t{decision}\t{}",
            by.as_str()
        )?;
    }
    Ok(())
}

// ── mbias.tsv / mbias_bounds.tsv ──────────────────────────────────────────────

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

/// Suggested mask lengths (from the CpG curve) per read/end.
fn write_mbias_bounds(
    w: &mut dyn Write,
    mbias: &MbiasAccumulator,
    detect: DetectParams,
    sample: &str,
) -> std::io::Result<()> {
    writeln!(w, "sample\tread\tend\tsuggested_mask\tplateau_fraction")?;
    for &role in &ReadRole::ALL {
        for &end in &ReadEnd::ALL {
            if !reported(role, end) || !mbias.has_data(role, end) {
                continue;
            }
            let k = detect_mask_length(mbias, role, end, detect);
            writeln!(
                w,
                "{sample}\t{}\t{}\t{k}\t{:.3}",
                role.label(),
                end.label(),
                detect.plateau_fraction,
            )?;
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

/// Resolve `sample`: explicit override wins, else comma-join unique `@RG SM:`.
fn resolve_sample(header: &Header, sample_override: Option<&str>) -> String {
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
    samples.into_iter().collect::<Vec<_>>().join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
    use crate::sieve::{DecidedBy, PerContextCounters, Stats};

    fn genome_with(counters: PerContextCounters, n_eval: u64, n_unconv: u64) -> Stats {
        let mut s = Stats::new(&[]);
        s.genome.n_templates = n_eval;
        s.genome.n_evaluated = n_eval;
        s.genome.n_unconverted = n_unconv;
        s.genome.counters = counters;
        s
    }

    fn render_summary(stats: &Stats, mbias: &MbiasAccumulator) -> String {
        let mut buf = Vec::new();
        write_summary(&mut buf, stats, mbias, "S1").unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn summary_has_read_column_and_all_row() {
        let stats = genome_with(PerContextCounters::default(), 0, 0);
        let mbias = MbiasAccumulator::new();
        let text = render_summary(&stats, &mbias);
        let mut lines = text.lines();
        let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(hdr[3], "read");
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(hdr.len(), row.len());
        assert_eq!(row[2], "genome");
        assert_eq!(row[3], "all");
    }

    #[test]
    fn per_read_rows_come_from_mbias() {
        let mut counters = PerContextCounters::default();
        counters.add_counts(Context::CpA, 1, 1000); // 1 unconverted of 1000
        let stats = genome_with(counters, 100, 1);
        let mut mbias = MbiasAccumulator::new();
        // R1: 2 unconverted of 10 CpA at cycle 0 → conv rate 0.8.
        for i in 0..10 {
            mbias.record(ReadRole::R1, ReadEnd::FivePrime, Context::CpA, 0, i < 2);
        }
        let text = render_summary(&stats, &mbias);
        let r1 = text.lines().find(|l| l.contains("\tR1\t")).expect("R1 row present");
        let cols: Vec<&str> = r1.split('\t').collect();
        assert_eq!(cols[9], "2", "CA_unconv");
        assert_eq!(cols[10], "10", "CA_total");
        assert_eq!(cols[4], "", "n_templates blank on per-read row");
        assert_eq!(cols[17], "0.800000", "conv_rate_CpA");
    }

    #[test]
    fn no_mbias_data_yields_only_all_rows() {
        let stats = genome_with(PerContextCounters::default(), 5, 0);
        let mbias = MbiasAccumulator::new();
        // header + genome/all only.
        assert_eq!(render_summary(&stats, &mbias).lines().count(), 2);
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
                (false, DecidedBy::ZeroSites)
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
        assert_eq!(lines[1], "S1\t0\t0\t\t5\tconverted\tzero_sites");
        assert_eq!(lines[2], "S1\t10\t3\t0.700000\t7\tunconverted\tcount");
        assert_eq!(lines[3], "S1\t50\t5\t0.900000\t2\tunconverted\tproportion");
    }
}
