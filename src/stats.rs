//! Multi-row stats TSV emitted when `--stats PATH` is set.
//!
//! One row per reporting scope: `genome` first (everything not on a control
//! contig), then one row per `--control-contig` in the order they were passed.
//! Per-context columns are broken out for **all four** contexts (CpA/CpC/CpT/CpG)
//! regardless of which subset drives the decision — so users can read, e.g.,
//! CpG retention on a methylated pUC19 control. Whole-run diagnostic counters
//! attach to the `genome` row only.
//!
//! A single `COLUMNS` array is the source of truth for both the header row and
//! every value row, so they cannot drift.

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

use anyhow::{Context as _, Result};
use noodles_sam::Header;
use noodles_sam::header::record::value::map::read_group::tag as rg_tag;

use crate::METHYLSIEVE_BUILD;
use crate::reference::Context;
use crate::sieve::{DecidedBy, ScopeStats, Stats};

/// Whole-run diagnostic counters, rendered only on the `genome` row.
#[derive(Debug, Clone, Copy)]
struct Diagnostics {
    chimeric_to_control_templates: u64,
    unmapped_templates: u64,
    zero_site_templates: u64,
    below_min_sites_templates: u64,
}

/// A single TSV row: scope stats plus shared context (sample, version) and
/// the optional whole-run diagnostics (present only on the genome row).
struct Row<'a> {
    sample: &'a str,
    scope: &'a ScopeStats,
    diagnostics: Option<Diagnostics>,
}

/// Per-context columns in report order: CpH first (A, C, T) then CpG.
const REPORT_ORDER: [Context; 4] = [Context::CpA, Context::CpC, Context::CpT, Context::CpG];

type ColumnFn = fn(&Row) -> String;
const COLUMNS: &[(&str, ColumnFn)] = &[
    ("sample", |r| r.sample.to_string()),
    ("methylsieve_version", |_| METHYLSIEVE_BUILD.to_string()),
    ("scope", |r| r.scope.name.clone()),
    ("n_templates", |r| r.scope.n_templates.to_string()),
    ("n_evaluated", |r| r.scope.n_evaluated.to_string()),
    ("n_unconverted", |r| r.scope.n_unconverted.to_string()),
    ("n_removed", |r| r.scope.n_removed.to_string()),
    ("frac_unconverted", |r| frac_unconverted(r.scope)),
    ("CA_unconv", |r| r.scope.counters.unconv[Context::CpA.index()].to_string()),
    ("CA_total", |r| r.scope.counters.total[Context::CpA.index()].to_string()),
    ("CC_unconv", |r| r.scope.counters.unconv[Context::CpC.index()].to_string()),
    ("CC_total", |r| r.scope.counters.total[Context::CpC.index()].to_string()),
    ("CT_unconv", |r| r.scope.counters.unconv[Context::CpT.index()].to_string()),
    ("CT_total", |r| r.scope.counters.total[Context::CpT.index()].to_string()),
    ("CG_unconv", |r| r.scope.counters.unconv[Context::CpG.index()].to_string()),
    ("CG_total", |r| r.scope.counters.total[Context::CpG.index()].to_string()),
    ("conv_rate_CpA", |r| conv_rate(r.scope, Context::CpA)),
    ("conv_rate_CpC", |r| conv_rate(r.scope, Context::CpC)),
    ("conv_rate_CpT", |r| conv_rate(r.scope, Context::CpT)),
    ("conv_rate_CpG", |r| conv_rate(r.scope, Context::CpG)),
    ("chimeric_to_control_templates", |r| {
        r.diagnostics.map(|d| d.chimeric_to_control_templates.to_string()).unwrap_or_default()
    }),
    ("unmapped_templates", |r| {
        r.diagnostics.map(|d| d.unmapped_templates.to_string()).unwrap_or_default()
    }),
    ("zero_site_templates", |r| {
        r.diagnostics.map(|d| d.zero_site_templates.to_string()).unwrap_or_default()
    }),
    ("below_min_sites_templates", |r| {
        r.diagnostics.map(|d| d.below_min_sites_templates.to_string()).unwrap_or_default()
    }),
];

/// Fraction of evaluated (evidence-bearing) templates decided unconverted.
/// Empty when no templates were evaluated.
fn frac_unconverted(scope: &ScopeStats) -> String {
    if scope.n_evaluated == 0 {
        String::new()
    } else {
        format!("{:.6}", scope.n_unconverted as f64 / scope.n_evaluated as f64)
    }
}

/// Conversion rate `1 - unconv/total` for `ctx`. Empty when no sites.
fn conv_rate(scope: &ScopeStats, ctx: Context) -> String {
    let total = scope.counters.total[ctx.index()];
    if total == 0 {
        String::new()
    } else {
        let unconv = scope.counters.unconv[ctx.index()];
        format!("{:.6}", 1.0 - unconv as f64 / total as f64)
    }
}

/// Render the full multi-row TSV (header + one row per scope) to `w`.
///
/// # Errors
/// Propagates write failures.
pub(crate) fn write_tsv<W: Write>(
    w: &mut W,
    stats: &Stats,
    header: &Header,
    sample_override: Option<&str>,
) -> Result<()> {
    // `REPORT_ORDER` is referenced for documentation/ordering intent; the
    // explicit per-context columns above hard-code the same order.
    debug_assert_eq!(REPORT_ORDER.len(), 4);

    let sample = resolve_sample(header, sample_override);
    let diagnostics = Diagnostics {
        chimeric_to_control_templates: stats.chimeric_to_control_templates,
        unmapped_templates: stats.unmapped_templates,
        zero_site_templates: stats.zero_site_templates,
        below_min_sites_templates: stats.below_min_sites_templates,
    };

    // Header row.
    for (i, (name, _)) in COLUMNS.iter().enumerate() {
        if i > 0 {
            w.write_all(b"\t")?;
        }
        w.write_all(name.as_bytes())?;
    }
    w.write_all(b"\n")?;

    // Genome row carries diagnostics; control rows do not.
    let rows = std::iter::once(Row {
        sample: &sample,
        scope: &stats.genome,
        diagnostics: Some(diagnostics),
    })
    .chain(stats.controls.iter().map(|c| Row {
        sample: &sample,
        scope: c,
        diagnostics: None,
    }));
    for row in rows {
        for (i, (_, render)) in COLUMNS.iter().enumerate() {
            if i > 0 {
                w.write_all(b"\t")?;
            }
            w.write_all(render(&row).as_bytes())?;
        }
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Open `path` (or stdout if `-`) and write the TSV.
///
/// # Errors
/// Propagates file/IO errors.
pub(crate) fn write_to_path(
    path: &Path,
    stats: &Stats,
    header: &Header,
    sample_override: Option<&str>,
) -> Result<()> {
    if path.to_string_lossy() == "-" {
        let mut stdout = std::io::stdout().lock();
        write_tsv(&mut stdout, stats, header, sample_override)
            .context("writing stats to stdout")?;
    } else {
        let mut f =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        write_tsv(&mut f, stats, header, sample_override)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

/// Render the conversion matrix — one row per observed `(checked, unconverted)`
/// decision cell — to `w`. `classify` replays the per-cell verdict (the
/// processor's [`crate::sieve::RecordProcessor::classify`]); it is pure in the
/// cell key, so the matrix never drifts from the live decision.
///
/// # Errors
/// Propagates write failures.
pub(crate) fn write_matrix_tsv<W, F>(
    w: &mut W,
    stats: &Stats,
    header: &Header,
    sample_override: Option<&str>,
    classify: F,
) -> Result<()>
where
    W: Write,
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    let sample = resolve_sample(header, sample_override);
    writeln!(
        w,
        "sample\tchecked_sites\tunconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
    )?;
    // BTreeMap key is (checked, unconverted), so rows come out sorted.
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

/// Open `path` (or stdout if `-`) and write the conversion matrix.
///
/// # Errors
/// Propagates file/IO errors.
pub(crate) fn write_matrix_to_path<F>(
    path: &Path,
    stats: &Stats,
    header: &Header,
    sample_override: Option<&str>,
    classify: F,
) -> Result<()>
where
    F: Fn(u64, u64) -> (bool, DecidedBy),
{
    if path.to_string_lossy() == "-" {
        let mut stdout = std::io::stdout().lock();
        write_matrix_tsv(&mut stdout, stats, header, sample_override, classify)
            .context("writing conversion matrix to stdout")?;
    } else {
        let mut f =
            std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
        write_matrix_tsv(&mut f, stats, header, sample_override, classify)
            .with_context(|| format!("writing {}", path.display()))?;
    }
    Ok(())
}

/// Resolve `sample`: explicit override wins, else comma-join the unique
/// `@RG SM:` tags, else empty.
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
    use crate::sieve::{DecidedBy, PerContextCounters, Stats};

    fn genome_with(counters: PerContextCounters, n_eval: u64, n_unconv: u64) -> Stats {
        let mut s = Stats::new(&[]);
        s.genome.n_templates = n_eval;
        s.genome.n_evaluated = n_eval;
        s.genome.n_unconverted = n_unconv;
        s.genome.counters = counters;
        s
    }

    #[test]
    fn header_and_values_have_equal_columns() {
        let stats = genome_with(PerContextCounters::default(), 0, 0);
        let header = Header::default();
        let mut buf = Vec::new();
        write_tsv(&mut buf, &stats, &header, Some("S1")).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let mut lines = text.lines();
        let hdr = lines.next().unwrap().split('\t').count();
        let val = lines.next().unwrap().split('\t').count();
        assert_eq!(hdr, COLUMNS.len());
        assert_eq!(hdr, val);
    }

    #[test]
    fn conv_rate_empty_when_no_sites() {
        let stats = genome_with(PerContextCounters::default(), 0, 0);
        let scope = &stats.genome;
        assert_eq!(conv_rate(scope, Context::CpG), "");
        assert_eq!(frac_unconverted(scope), "");
    }

    #[test]
    fn conv_rate_computed_from_counts() {
        let mut counters = PerContextCounters::default();
        // 1 unconverted of 1000 CpA → conv rate 0.999.
        counters.total[Context::CpA.index()] = 1000;
        counters.unconv[Context::CpA.index()] = 1;
        let stats = genome_with(counters, 100, 1);
        let scope = &stats.genome;
        assert_eq!(conv_rate(scope, Context::CpA), "0.999000");
        assert_eq!(frac_unconverted(scope), "0.010000");
    }

    #[test]
    fn control_rows_follow_genome_and_omit_diagnostics() {
        let mut stats = Stats::new(&["phage_lambda".to_string()]);
        stats.genome.n_templates = 10;
        stats.genome.n_evaluated = 10;
        stats.unmapped_templates = 3;
        stats.controls[0].n_templates = 5;
        let header = Header::default();
        let mut buf = Vec::new();
        write_tsv(&mut buf, &stats, &header, Some("S1")).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3); // header + genome + 1 control
        assert!(lines[1].contains("genome"));
        assert!(lines[2].contains("phage_lambda"));
        // The last four columns are the diagnostics (chimeric, unmapped,
        // zero_site, below_min_sites): populated on the genome row, blank on the
        // control row.
        let genome_cols: Vec<&str> = lines[1].split('\t').collect();
        let control_cols: Vec<&str> = lines[2].split('\t').collect();
        let n = COLUMNS.len();
        assert_eq!(&genome_cols[n - 4..], &["0", "3", "0", "0"]);
        assert_eq!(&control_cols[n - 4..], &["", "", "", ""]);
    }

    #[test]
    fn conversion_matrix_renders_sorted_cells_with_verdicts() {
        let mut stats = Stats::new(&[]);
        // (checked_sites, unconverted_sites) → n_templates, inserted out of order.
        stats.conversion_matrix.insert((50, 5), 2);
        stats.conversion_matrix.insert((0, 0), 5);
        stats.conversion_matrix.insert((10, 3), 7);
        // Stand-in for RecordProcessor::classify: count arm below 40 sites,
        // proportion arm at/above, zero sites → converted.
        let classify = |unconv: u64, checked: u64| -> (bool, DecidedBy) {
            if checked == 0 {
                (false, DecidedBy::ZeroSites)
            } else if checked >= 40 {
                (unconv as f64 / checked as f64 > 0.05, DecidedBy::Proportion)
            } else {
                (unconv >= 3, DecidedBy::Count)
            }
        };
        let header = Header::default();
        let mut buf = Vec::new();
        write_matrix_tsv(&mut buf, &stats, &header, Some("S1"), classify).unwrap();
        let lines: Vec<String> =
            String::from_utf8(buf).unwrap().lines().map(str::to_string).collect();
        assert_eq!(
            lines[0],
            "sample\tchecked_sites\tunconverted_sites\tconversion_rate\tn_templates\tdecision\tdecided_by"
        );
        // Sorted by (checked, unconverted); conv_rate blank when no sites.
        assert_eq!(lines[1], "S1\t0\t0\t\t5\tconverted\tzero_sites");
        assert_eq!(lines[2], "S1\t10\t3\t0.700000\t7\tunconverted\tcount");
        assert_eq!(lines[3], "S1\t50\t5\t0.900000\t2\tunconverted\tproportion");
        assert_eq!(lines.len(), 4);
    }
}
