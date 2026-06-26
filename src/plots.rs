//! PDF plots for `--metrics-prefix`, rendered with kuva from the same in-memory
//! data as the metric TSVs (see [`crate::metrics`]):
//!
//! - `PREFIX.mbias.pdf` — per-read-cycle cytosine-retention curves by context
//!   (one panel per read), with the applied 5' mask shaded.
//! - `PREFIX.conversion-matrix.pdf` — a hexbin of per-template observed-vs-
//!   converted sites (log-scaled template density) with the converted/unconverted
//!   decision boundary overlaid.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context as _, Result, anyhow};
use kuva::backend::pdf::PdfBackend;
use kuva::plot::hexbin::{HexbinPlot, ZReduce};
use kuva::plot::legend::LegendPosition;
use kuva::plot::{ColorMap, LinePlot};
use kuva::render::annotations::{ReferenceLine, ShadedRegion};
use kuva::render::figure::Figure;
use kuva::render::layout::{Layout, TickFormat};
use kuva::render::plots::Plot;
use kuva::render::render::{Scene, render_multiple};

use crate::mask::MaskPlan;
use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
use crate::reference::Context;
use crate::sieve::DecidedBy;

// Fulcrum Genomics brand colors (source: riker `plotting.rs`).
const FG_BLUE: &str = "#26a8e0";
const FG_TEAL: &str = "#2fae99";
const FG_FOREST: &str = "#269e2a";
const FG_RED: &str = "#e04040";
const FG_GRAY: &str = "#5c7682";
/// Bright cyan for the decision boundary — absent from the Inferno colormap, so
/// it stays legible over both the dark and bright hexes.
const BOUNDARY: &str = "#00e5ff";
/// Purple for the M-bias detection threshold — distinct from the four
/// context colors (A/C/G/T) and the gray mask marks.
const FG_PURPLE: &str = "#7c3aed";

/// Per-context line color: CpG (the methylated context) salient red, CpH cool.
fn ctx_color(ctx: Context) -> &'static str {
    match ctx {
        Context::CpA => FG_BLUE,
        Context::CpC => FG_TEAL,
        Context::CpT => FG_FOREST,
        Context::CpG => FG_RED,
    }
}

/// Render a built scene to a PDF file.
fn write_scene_pdf(path: &Path, scene: &Scene) -> Result<()> {
    let bytes =
        PdfBackend.render_scene(scene).map_err(|e| anyhow!("rendering {}: {e}", path.display()))?;
    std::fs::write(path, bytes).with_context(|| format!("writing {}", path.display()))
}

// ── M-bias ──────────────────────────────────────────────────────────────────

/// One read's panel: a retention-rate-vs-cycle line per context (5' end).
fn mbias_panel(mbias: &MbiasAccumulator, role: ReadRole, legend: bool) -> Vec<Plot> {
    Context::ALL
        .iter()
        .filter_map(|&ctx| {
            let pts: Vec<(f64, f64)> = mbias
                .cycles(role, ReadEnd::FivePrime, ctx)
                .iter()
                .enumerate()
                .filter_map(|(i, cc)| cc.frac().map(|f| ((i + 1) as f64, f)))
                .collect();
            if pts.is_empty() {
                return None;
            }
            let mut lp =
                LinePlot::new().with_data(pts).with_color(ctx_color(ctx)).with_stroke_width(1.8);
            if legend {
                lp = lp.with_legend(ctx.label());
            }
            Some(Plot::Line(lp))
        })
        .collect()
}

/// One read's panel layout, shading the applied 5' mask when masking ran and
/// drawing the `fraction × plateau` detection threshold as a horizontal cut line.
fn mbias_layout(
    title: &str,
    y_label: bool,
    max_cycle: f64,
    mask: Option<usize>,
    threshold: Option<f64>,
    plateau_fraction: f64,
) -> Layout {
    let mut l = Layout::new((0.0, max_cycle), (0.0, 1.0))
        .with_title(title)
        .with_x_label("Position in Read")
        // Stop the x-axis at the read length rather than a major tick wider.
        .with_x_axis_min(0.0)
        .with_x_axis_max(max_cycle)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        .with_show_grid(true)
        // Smaller tick + reference-line labels (the "mask N" / threshold labels).
        .with_tick_size(9);
    if y_label {
        l = l.with_y_label("Cytosine Retention Rate");
    }
    if let Some(k) = mask.filter(|k| *k > 0) {
        let mut region = ShadedRegion::vertical(0.5, k as f64 + 0.5);
        region.color = FG_GRAY.to_string();
        region.opacity = 0.10;
        l = l.with_shaded_region(region).with_reference_line(
            ReferenceLine::vertical(k as f64 + 0.5)
                .with_color(FG_GRAY)
                .with_stroke_width(1.0)
                .with_dasharray("4,3")
                .with_label(format!("mask {k}")),
        );
    }
    // The horizontal detection threshold: masking cuts just before the smoothed
    // CpG curve first holds above this level. A thin solid purple line; its label
    // is right-aligned to the read-length edge by the renderer (End anchor).
    if let Some(thr) = threshold.filter(|t| t.is_finite() && *t > 0.0) {
        l = l.with_reference_line(
            ReferenceLine::horizontal(thr)
                .with_color(FG_PURPLE)
                .with_stroke_width(0.8)
                .with_dasharray("none")
                .with_label(format!("{plateau_fraction:.2} × plateau")),
        );
    }
    l
}

/// Write the M-bias curves to `path`. No-op (no file) when no read class has
/// data. `mask_plan` shades the applied 5' mask per read; `None` when masking
/// was not run.
pub(crate) fn write_mbias_pdf(
    path: &Path,
    mbias: &MbiasAccumulator,
    mask_plan: Option<&MaskPlan>,
    sample: &str,
) -> Result<()> {
    let roles: Vec<ReadRole> =
        ReadRole::ALL.iter().copied().filter(|&r| mbias.has_data(r, ReadEnd::FivePrime)).collect();
    if roles.is_empty() {
        return Ok(());
    }
    let max_cycle = roles
        .iter()
        .flat_map(|&r| {
            Context::ALL.iter().map(move |&c| mbias.cycles(r, ReadEnd::FivePrime, c).len())
        })
        .max()
        .unwrap_or(1) as f64;
    let role_title = |r: ReadRole| match r {
        ReadRole::R1 => "Read 1",
        ReadRole::R2 => "Read 2",
        ReadRole::Se => "Single-End",
    };

    let frac = mask_plan.map_or(0.0, MaskPlan::plateau_fraction);
    let mut panels = Vec::with_capacity(roles.len());
    let mut layouts = Vec::with_capacity(roles.len());
    for (i, &role) in roles.iter().enumerate() {
        let first = i == 0;
        panels.push(mbias_panel(mbias, role, first));
        layouts.push(mbias_layout(
            role_title(role),
            first,
            max_cycle,
            mask_plan.map(|p| p.role_5p(role)),
            mask_plan.and_then(|p| p.threshold_5p(role)),
            frac,
        ));
    }
    let scene = Figure::new(1, roles.len())
        .with_plots(panels)
        .with_layouts(layouts)
        .with_shared_legend()
        .with_title(format!("M-bias — {sample}"))
        .with_cell_size(480.0, 380.0)
        .render();
    write_scene_pdf(path, &scene)
}

// ── Conversion matrix ─────────────────────────────────────────────────────────

/// Write the conversion-matrix hexbin to `path`. `matrix` is keyed by
/// `(checked, unconverted)` → template count over the decision contexts;
/// `classify` replays the per-cell verdict for the boundary; `contexts` is the
/// short context label (e.g. `"CpH"`). No-op when the matrix is empty.
pub(crate) fn write_matrix_pdf(
    path: &Path,
    matrix: &BTreeMap<(u64, u64), u64>,
    classify: impl Fn(u64, u64) -> (bool, DecidedBy),
    sample: &str,
    contexts: &str,
) -> Result<()> {
    if matrix.is_empty() {
        return Ok(());
    }
    // Template totals for the title summary: how many templates the decision
    // flags as unconverted, out of all that hit the matrix.
    let total: u64 = matrix.values().sum();
    let unconverted: u64 =
        matrix.iter().filter_map(|(&(c, u), &n)| classify(u, c).0.then_some(n)).sum();

    // Clamp the axes to the 99.9th percentile of observed sites. A handful of
    // supplementary-inflated outliers (split/chimeric reads mapping to several
    // loci, each adding sites) otherwise stretch the range and squash the bulk
    // of the data into a corner.
    let cap = percentile_checked(matrix, total, 0.999);
    let span = cap as f64 + 1.0;

    let (mut xs, mut ys, mut zs) = (Vec::new(), Vec::new(), Vec::new());
    // Lowest converted-count still called "converted" per observed count → the
    // decision boundary drawn over the hexbin.
    let mut boundary: BTreeMap<u64, u64> = BTreeMap::new();
    for (&(checked, unconv), &n) in matrix {
        if checked > cap {
            continue; // beyond the clamped range (rare outliers)
        }
        let converted_sites = checked - unconv;
        xs.push(checked as f64);
        ys.push(converted_sites as f64);
        zs.push(n as f64);
        if !classify(unconv, checked).0 {
            boundary
                .entry(checked)
                .and_modify(|m| *m = (*m).min(converted_sites))
                .or_insert(converted_sites);
        }
    }

    let hex = HexbinPlot::new()
        .with_data(xs, ys)
        .with_z(zs, ZReduce::Sum)
        .with_log_color(true)
        .with_color_map(ColorMap::Inferno)
        .with_n_bins((cap as usize).max(1))
        .with_x_range(0.0, span)
        .with_y_range(0.0, span)
        .with_min_count(1)
        .with_stroke("#ffffff")
        .with_stroke_width(0.4)
        .with_colorbar(true)
        .with_colorbar_label("Templates (log)");
    let bpts: Vec<(f64, f64)> =
        boundary.iter().map(|(&x, &c)| (x as f64, c as f64 - 0.5)).collect();
    let bound = LinePlot::new()
        .with_data(bpts)
        .with_color(BOUNDARY)
        .with_stroke_width(0.8)
        .with_step()
        .with_legend("converted / unconverted boundary");

    // Template-split summary as a second title line (kuva has no subtitle; a `\n`
    // in the title word-wraps into another line). The blank middle line widens the
    // gap — kuva spaces title lines at exactly 1× font size, which otherwise reads
    // as jammed. Both lines render at title size pending an upstream subtitle.
    let pct = if total > 0 { 100.0 * unconverted as f64 / total as f64 } else { 0.0 };
    let title = format!(
        "{contexts} Conversion in Templates from {sample}\n\n{} Templates · {} ({pct:.1}%) Flagged as Unconverted",
        count_label(total as f64),
        count_label(unconverted as f64),
    );

    let plots = vec![Plot::Hexbin(hex), Plot::Line(bound)];
    // auto_from_plots so the hexbin colorbar is detected; force a square range.
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_title_wrap(80) // splits on the embedded newline; lines stay under 80
        .with_x_label(format!("{contexts} Sites Observed"))
        .with_y_label("Converted Sites")
        .with_x_axis_min(0.0)
        .with_x_axis_max(span)
        .with_y_axis_min(0.0)
        .with_y_axis_max(span)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        // SI-compact colorbar ticks ("100k") so 6-digit counts don't clip/overlap.
        .with_colorbar_tick_format(TickFormat::Custom(Arc::new(|v: f64| si_compact(v))))
        // Frameless legend laid over the grid (no box/outline).
        .with_legend_position(LegendPosition::InsideTopLeft)
        .with_legend_box(false);
    let scene = render_multiple(plots, layout);
    write_scene_pdf(path, &scene)
}

/// SI-compact number for colorbar tick labels: `100000 → "100k"`, `1.2e6 → "1.2M"`,
/// smaller values rounded to an integer.
fn si_compact(v: f64) -> String {
    let a = v.abs();
    if a >= 1e6 {
        format!("{:.1}M", v / 1e6)
    } else if a >= 1e3 {
        format!("{:.0}k", v / 1e3)
    } else {
        format!("{v:.0}")
    }
}

/// Compact count for the title summary: `3.0e6 → "3.0m"`, `90000 → "90K"`.
fn count_label(v: f64) -> String {
    let a = v.abs();
    if a >= 1e6 {
        format!("{:.1}m", v / 1e6)
    } else if a >= 1e3 {
        format!("{:.0}K", v / 1e3)
    } else {
        format!("{v:.0}")
    }
}

/// Smallest observed-site count `C` (rounded up to a multiple of 10) such that
/// templates with ≤ `C` sites cover `frac` of all templates — used to clamp the
/// axis past the thin tail of supplementary-inflated outliers.
fn percentile_checked(matrix: &BTreeMap<(u64, u64), u64>, total: u64, frac: f64) -> u64 {
    if total == 0 {
        return 10;
    }
    let target = (total as f64 * frac).ceil() as u64;
    let mut by_checked: BTreeMap<u64, u64> = BTreeMap::new();
    for (&(c, _), &n) in matrix {
        *by_checked.entry(c).or_default() += n;
    }
    let mut acc = 0u64;
    for (&c, &n) in &by_checked {
        acc += n;
        if acc >= target {
            return c.div_ceil(10) * 10;
        }
    }
    by_checked.keys().last().copied().unwrap_or(1).div_ceil(10) * 10
}
