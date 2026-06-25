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

use anyhow::{Context as _, Result, anyhow};
use kuva::backend::pdf::PdfBackend;
use kuva::plot::hexbin::{HexbinPlot, ZReduce};
use kuva::plot::legend::LegendPosition;
use kuva::plot::{ColorMap, LinePlot};
use kuva::render::annotations::{ReferenceLine, ShadedRegion};
use kuva::render::figure::Figure;
use kuva::render::layout::Layout;
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

/// One read's panel layout, shading the applied 5' mask when masking ran.
fn mbias_layout(title: &str, y_label: bool, max_cycle: f64, mask: Option<usize>) -> Layout {
    let mut l = Layout::new((0.0, max_cycle), (0.0, 1.0))
        .with_title(title)
        .with_x_label("Position in Read")
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        .with_show_grid(true);
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
    let max_checked = matrix.keys().map(|&(c, _)| c).max().unwrap_or(0);
    let span = max_checked as f64 + 1.0;

    let (mut xs, mut ys, mut zs) = (Vec::new(), Vec::new(), Vec::new());
    // Lowest converted-count still called "converted" per observed count → the
    // decision boundary (a position masked in one mate is masked in the other...).
    let mut boundary: BTreeMap<u64, u64> = BTreeMap::new();
    for (&(checked, unconv), &n) in matrix {
        let converted = checked - unconv;
        xs.push(checked as f64);
        ys.push(converted as f64);
        zs.push(n as f64);
        if !classify(unconv, checked).0 {
            boundary.entry(checked).and_modify(|m| *m = (*m).min(converted)).or_insert(converted);
        }
    }

    let hex = HexbinPlot::new()
        .with_data(xs, ys)
        .with_z(zs, ZReduce::Sum)
        .with_log_color(true)
        .with_color_map(ColorMap::Inferno)
        .with_n_bins((max_checked as usize).max(1))
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

    let plots = vec![Plot::Hexbin(hex), Plot::Line(bound)];
    // auto_from_plots so the hexbin colorbar is detected; force a square range.
    let layout = Layout::auto_from_plots(&plots)
        .with_title(format!("{contexts} Conversion in Templates from {sample}"))
        .with_x_label(format!("{contexts} Sites Observed"))
        .with_y_label("Converted Sites")
        .with_x_axis_min(0.0)
        .with_x_axis_max(span)
        .with_y_axis_min(0.0)
        .with_y_axis_max(span)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        .with_legend_position(LegendPosition::InsideTopLeft)
        .with_legend_width(210.0)
        .with_legend_height(30.0);
    let scene = render_multiple(plots, layout);
    write_scene_pdf(path, &scene)
}
