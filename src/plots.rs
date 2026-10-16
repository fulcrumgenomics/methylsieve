//! PDF plots for `--metrics-prefix`, rendered with kuva from the same in-memory
//! data as the metric TSVs (see [`crate::metrics`]):
//!
//! - `PREFIX.mbias.pdf` — per-read-cycle apparent-methylation curves by context
//!   (one panel per read end), read through the library chemistry, with the
//!   applied mask shaded — the 5' mask for every read, plus the single-end 3'
//!   mask on its own panel.
//! - `PREFIX.conversion-matrix.pdf` — one page per scope (the genome, then each
//!   control): a hexbin of per-template observed-vs-converted sites (log-scaled
//!   template density). Chemistries that make the unconverted call overlay the
//!   converted/unconverted boundary; the others draw the `x = y` diagonal of
//!   complete conversion instead, so the sag below it reads at a glance.

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

use crate::chemistry::Chemistry;
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
/// Bright green for the complete-conversion diagonal — distinct from the cyan
/// boundary and, like it, absent from the Inferno colormap.
const DIAGONAL: &str = "#76ff03";

/// Per-context line color: CpG (the methylated context) salient red, CpH cool.
fn ctx_color(ctx: Context) -> &'static str {
    match ctx {
        Context::CpA => FG_BLUE,
        Context::CpC => FG_TEAL,
        Context::CpT => FG_FOREST,
        Context::CpG => FG_RED,
    }
}

/// Render built scenes to a PDF file, one page per scene.
fn write_scenes_pdf(path: &Path, scenes: &[Scene]) -> Result<()> {
    let bytes = PdfBackend::new()
        .render_scenes(scenes)
        .map_err(|e| anyhow!("rendering {}: {e}", path.display()))?;
    std::fs::write(path, bytes).with_context(|| format!("writing {}", path.display()))
}

// ── M-bias ──────────────────────────────────────────────────────────────────

/// One read's panel: a methylation-vs-cycle line per context, for the given
/// read `end` (5' cycles counted from the read start, 3' from the read end).
fn mbias_panel(mbias: &MbiasAccumulator, role: ReadRole, end: ReadEnd, legend: bool) -> Vec<Plot> {
    let readout = mbias.readout();
    Context::ALL
        .iter()
        .filter_map(|&ctx| {
            let pts: Vec<(f64, f64)> = mbias
                .cycles(role, end, ctx)
                .iter()
                .enumerate()
                .filter_map(|(i, cc)| cc.frac_methylated(readout).map(|f| ((i + 1) as f64, f)))
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
    x_label: &str,
    y_label: bool,
    max_cycle: f64,
    mask: Option<usize>,
    threshold: Option<f64>,
    plateau_fraction: f64,
) -> Layout {
    let mut l = Layout::new((0.0, max_cycle), (0.0, 1.0))
        .with_title(title)
        .with_x_label(x_label)
        // Stop the x-axis at the read length rather than a major tick wider.
        .with_x_axis_min(0.0)
        .with_x_axis_max(max_cycle)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        .with_show_grid(true)
        // Smaller tick + reference-line labels (the "mask N" / threshold labels).
        .with_tick_size(9);
    if y_label {
        l = l.with_y_label("Apparent Methylation");
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
/// was not run. `chemistry` is named in the title.
pub(crate) fn write_mbias_pdf(
    path: &Path,
    mbias: &MbiasAccumulator,
    mask_plan: Option<&MaskPlan>,
    sample: &str,
    chemistry: Chemistry,
) -> Result<()> {
    let role_title = |r: ReadRole| match r {
        ReadRole::R1 => "Read 1",
        ReadRole::R2 => "Read 2",
        ReadRole::Se => "Single-End",
    };

    // One panel per (role, end) with data: the 5' curve for every role, plus the
    // single-end 3' curve, which has its own learned mask (`k_se_3p`) that the 5'
    // panels never show. The 3' panel plots position-from-3'-end, so its mask
    // shades from the left exactly like the 5' panels.
    struct PanelSpec {
        role: ReadRole,
        end: ReadEnd,
        title: String,
        x_label: &'static str,
        mask: Option<usize>,
        threshold: Option<f64>,
    }
    let mut specs: Vec<PanelSpec> = Vec::new();
    for &role in &ReadRole::ALL {
        if mbias.has_data(role, ReadEnd::FivePrime) {
            specs.push(PanelSpec {
                role,
                end: ReadEnd::FivePrime,
                title: role_title(role).to_string(),
                x_label: "Position in Read",
                mask: mask_plan.map(|p| p.role_5p(role)),
                threshold: mask_plan.and_then(|p| p.threshold_5p(role)),
            });
        }
    }
    if mbias.has_data(ReadRole::Se, ReadEnd::ThreePrime) {
        specs.push(PanelSpec {
            role: ReadRole::Se,
            end: ReadEnd::ThreePrime,
            title: "Single-End 3'".to_string(),
            x_label: "Position from 3' End",
            mask: mask_plan.map(|p| p.k_se_3p()),
            threshold: None, // no 3' plateau is stored; the applied mask is shaded
        });
    }
    if specs.is_empty() {
        return Ok(());
    }

    let max_cycle = specs
        .iter()
        .flat_map(|s| Context::ALL.iter().map(move |&c| mbias.cycles(s.role, s.end, c).len()))
        .max()
        .unwrap_or(1) as f64;

    let frac = mask_plan.map_or(0.0, MaskPlan::plateau_fraction);
    let mut panels = Vec::with_capacity(specs.len());
    let mut layouts = Vec::with_capacity(specs.len());
    for (i, s) in specs.iter().enumerate() {
        let first = i == 0;
        panels.push(mbias_panel(mbias, s.role, s.end, first));
        layouts.push(mbias_layout(
            &s.title,
            s.x_label,
            first,
            max_cycle,
            s.mask,
            s.threshold,
            frac,
        ));
    }
    let scene = Figure::new(1, specs.len())
        .with_plots(panels)
        .with_layouts(layouts)
        .with_shared_legend()
        .with_title(format!("M-bias — {sample} ({})", chemistry.name()))
        .with_cell_size(480.0, 380.0)
        .render();
    write_scenes_pdf(path, std::slice::from_ref(&scene))
}

// ── Conversion matrix ─────────────────────────────────────────────────────────

/// One page of the conversion-matrix PDF: a scope's per-template histogram.
pub(crate) struct MatrixPage<'a> {
    /// Scope name: `"genome"` or a control contig.
    pub(crate) scope: &'a str,
    /// Whether the scope is a spike-in control (named in the page title).
    pub(crate) is_control: bool,
    /// Short label of the contexts the histogram counts (e.g. `"CpH"`).
    pub(crate) contexts: &'a str,
    /// `(checked, unconverted)` → template count.
    pub(crate) cells: &'a BTreeMap<(u64, u64), u64>,
    /// Draw the converted/unconverted boundary (and count flagged templates);
    /// otherwise draw the complete-conversion diagonal.
    pub(crate) boundary: bool,
}

/// Write the conversion matrices to `path`, one page per non-empty page in
/// `pages`. `classify` replays the per-cell verdict for pages with a boundary.
/// No-op (no file) when every page is empty.
pub(crate) fn write_matrix_pdf(
    path: &Path,
    pages: &[MatrixPage<'_>],
    classify: impl Fn(u64, u64) -> (bool, DecidedBy),
    sample: &str,
) -> Result<()> {
    let scenes: Vec<Scene> = pages
        .iter()
        .filter(|p| !p.cells.is_empty())
        .map(|p| matrix_scene(p, &classify, sample))
        .collect();
    if scenes.is_empty() {
        return Ok(());
    }
    write_scenes_pdf(path, &scenes)
}

/// Sub-points per axis that each whole-number matrix cell is spread over before
/// hexbinning (so `SPREAD²` per cell). See [`spread_cell`].
const SPREAD: usize = 4;

/// A cell's unit square overlaps at most about four hexes. A cell with fewer
/// templates than that is drawn as one point: spreading it would only scatter
/// fractions of a template across several hexes, turning a lone template into a
/// cluster of faint ones.
const MIN_TEMPLATES_TO_SPREAD: u64 = 4;

/// Render one scope's hexbin page.
fn matrix_scene(
    page: &MatrixPage<'_>,
    classify: &impl Fn(u64, u64) -> (bool, DecidedBy),
    sample: &str,
) -> Scene {
    let matrix = page.cells;
    let contexts = page.contexts;
    // Template totals for the title summary: how many templates fall past the
    // boundary (the decision flags them on the genome), out of all that hit the
    // matrix.
    let total: u64 = matrix.values().sum();
    let unconverted: u64 = if page.boundary {
        matrix.iter().filter_map(|(&(c, u), &n)| classify(u, c).0.then_some(n)).sum()
    } else {
        0
    };

    // Clamp the axes to the 99.9th percentile of observed sites. A handful of
    // supplementary-inflated outliers (split/chimeric reads mapping to several
    // loci, each adding sites) otherwise stretch the range and squash the bulk
    // of the data into a corner.
    let cap = percentile_checked(matrix, total, 0.999);
    let span = cap as f64 + 1.0;

    let (mut xs, mut ys, mut zs) = (Vec::new(), Vec::new(), Vec::new());
    for (&(checked, unconv), &n) in matrix {
        if checked > cap {
            continue; // beyond the clamped range (rare outliers)
        }
        for (x, y, z) in spread_cell(checked, checked - unconv, n) {
            xs.push(x);
            ys.push(y);
            zs.push(z);
        }
    }
    let boundary = if page.boundary { decision_boundary(cap, classify) } else { BTreeMap::new() };

    let hex = HexbinPlot::new()
        .with_data(xs, ys)
        .with_z(zs, ZReduce::Sum)
        .with_log_color(true)
        .with_color_map(ColorMap::Inferno)
        // One bin per checked-site value up to a cap, so a heavy-tailed dataset
        // (thousands of sites on long templates) doesn't render a multi-thousand-
        // bin hexbin for a title/summary plot.
        .with_n_bins((cap as usize).clamp(1, 300))
        .with_x_range(-1.0, span)
        .with_y_range(-1.0, span)
        .with_min_count(1)
        .with_stroke("#ffffff")
        .with_stroke_width(0.4)
        .with_colorbar(true)
        .with_colorbar_label("Templates (log)");
    let overlay = if page.boundary {
        let bpts: Vec<(f64, f64)> =
            boundary.iter().map(|(&x, &c)| (x as f64, c as f64 - 0.5)).collect();
        LinePlot::new()
            .with_data(bpts)
            .with_color(BOUNDARY)
            .with_stroke_width(0.8)
            .with_step()
            .with_legend("converted / unconverted boundary")
    } else {
        LinePlot::new()
            .with_data(vec![(0.0, 0.0), (span, span)])
            .with_color(DIAGONAL)
            .with_stroke_width(0.8)
            .with_dashed()
            .with_legend("complete conversion (x = y)")
    };

    // Template-split summary as a second title line (kuva has no subtitle; a `\n`
    // in the title word-wraps into another line). The blank middle line widens the
    // gap — kuva spaces title lines at exactly 1× font size, which otherwise reads
    // as jammed. Both lines render at title size pending an upstream subtitle.
    let what =
        if page.is_control { format!("{} Templates", page.scope) } else { "Templates".into() };
    let summary = if page.boundary {
        let pct = if total > 0 { 100.0 * unconverted as f64 / total as f64 } else { 0.0 };
        // Control templates are never flagged, only counted against the boundary.
        let verdict = if page.is_control { "Past the Boundary" } else { "Flagged as Unconverted" };
        format!(
            "{} Templates · {} ({pct:.1}%) {verdict}",
            si_compact(total as f64),
            si_compact(unconverted as f64),
        )
    } else {
        format!("{} Templates", si_compact(total as f64))
    };
    let title = format!("{contexts} Conversion in {what} from {sample}\n\n{summary}");

    let plots = vec![Plot::Hexbin(hex), Plot::Line(overlay)];
    // auto_from_plots so the hexbin colorbar is detected; force a square range.
    let layout = Layout::auto_from_plots(&plots)
        .with_title(title)
        .with_title_wrap(80) // splits on the embedded newline; lines stay under 80
        .with_x_label(format!("{contexts} Sites Observed"))
        .with_y_label("Converted Sites")
        // Start both axes at -1 so the hexes centred on 0 (templates with no
        // converted sites, the failures) draw whole instead of being cut in half
        // by the axis. Major ticks fall on multiples of the step, so labels still
        // start at 0.
        .with_x_axis_min(-1.0)
        .with_x_axis_max(span)
        .with_y_axis_min(-1.0)
        .with_y_axis_max(span)
        .with_minor_ticks(5)
        .with_show_minor_grid(true)
        // SI-compact colorbar ticks ("100k") so 6-digit counts don't clip/overlap.
        .with_colorbar_tick_format(TickFormat::Custom(Arc::new(|v: f64| si_compact(v))))
        // Frameless legend laid over the grid (no box/outline).
        .with_legend_position(LegendPosition::InsideTopLeft)
        .with_legend_box(false);
    render_multiple(plots, layout)
}

/// The converted/unconverted boundary drawn over the hexbin: for every site count
/// from 0 to `cap`, the fewest converted sites `classify` still calls converted.
/// Traced from the decision rather than the observed cells, so it runs unbroken
/// across site counts no template has, and shows even when every template is past it.
fn decision_boundary(
    cap: u64,
    classify: &impl Fn(u64, u64) -> (bool, DecidedBy),
) -> BTreeMap<u64, u64> {
    (0..=cap)
        .filter_map(|checked| {
            let lowest = (0..=checked).find(|&conv| !classify(checked - conv, checked).0)?;
            Some((checked, lowest))
        })
        .collect()
}

/// The weighted points `(x, y, templates)` that represent the matrix cell of `n`
/// templates at `(checked, converted)` in the hexbin.
///
/// The counts sit on a whole-number lattice and the hex grid can't line up with
/// it, so one point per cell leaves some hexes with no cell (white holes in dense
/// regions) and sums neighbouring cells into others. Spreading the cell's
/// templates evenly over its unit square as a `SPREAD`×`SPREAD` grid of
/// sub-points instead gives each hex templates in proportion to the area it
/// covers. Cells below [`MIN_TEMPLATES_TO_SPREAD`] stay a single point. Nothing
/// lies below zero on either axis, so a cell on a zero edge is not spread across
/// it: the row (or column) at zero stays one row of hexes.
fn spread_cell(checked: u64, converted: u64, n: u64) -> Vec<(f64, f64, f64)> {
    let (x, y) = (checked as f64, converted as f64);
    if n < MIN_TEMPLATES_TO_SPREAD {
        return vec![(x, y, n as f64)];
    }
    // Sub-point offsets across the unit cell: -0.375, -0.125, 0.125, 0.375 for four.
    let offsets = |at_zero: bool| -> Vec<f64> {
        if at_zero {
            vec![0.0]
        } else {
            (0..SPREAD).map(|i| (i as f64 + 0.5) / SPREAD as f64 - 0.5).collect()
        }
    };
    let (dxs, dys) = (offsets(checked == 0), offsets(converted == 0));
    let share = n as f64 / (dxs.len() * dys.len()) as f64;
    dxs.iter().flat_map(|&dx| dys.iter().map(move |&dy| (x + dx, y + dy, share))).collect()
}

/// SI-compact number for colorbar tick labels and the title summary:
/// `100000 → "100k"`, `1.2e6 → "1.2M"`, smaller values rounded to an integer.
/// Uses conventional SI casing (`k` = 10³, `M` = 10⁶).
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

#[cfg(test)]
mod tests {
    use super::*;

    /// More than `max_unconverted` unconverted sites flags a template.
    fn count_rule(max_unconverted: u64) -> impl Fn(u64, u64) -> (bool, DecidedBy) {
        move |unconv, _checked| (unconv > max_unconverted, DecidedBy::Count)
    }

    #[test]
    fn decision_boundary_covers_every_site_count() {
        let boundary = decision_boundary(10, &count_rule(3));
        assert_eq!(boundary.keys().copied().collect::<Vec<_>>(), (0..=10).collect::<Vec<_>>());
    }

    #[test]
    fn decision_boundary_is_the_fewest_converted_sites_still_converted() {
        let boundary = decision_boundary(10, &count_rule(3));
        assert_eq!(boundary[&2], 0, "too few sites to flag at all");
        assert_eq!(boundary[&10], 7, "at most 3 of 10 unconverted");
    }

    #[test]
    fn spread_cell_does_not_spread_below_zero_converted() {
        let points = spread_cell(20, 0, 1000);
        assert!(points.iter().all(|&(_, y, _)| y == 0.0));
        let total: f64 = points.iter().map(|&(_, _, z)| z).sum();
        assert!((total - 1000.0).abs() < 1e-9, "spread total {total}");
    }

    #[test]
    fn spread_cell_keeps_zero_site_templates_at_the_origin() {
        assert_eq!(spread_cell(0, 0, 1000), vec![(0.0, 0.0, 1000.0)]);
    }

    #[test]
    fn spread_cell_keeps_every_template() {
        let total: f64 = spread_cell(20, 7, 1000).iter().map(|&(_, _, z)| z).sum();
        assert!((total - 1000.0).abs() < 1e-9, "spread total {total}");
    }

    #[test]
    fn spread_cell_stays_within_its_unit_square() {
        for (x, y, _) in spread_cell(20, 7, 1000) {
            assert!((19.5..20.5).contains(&x) && (6.5..7.5).contains(&y), "({x}, {y})");
        }
    }

    #[test]
    fn spread_cell_spreads_cells_at_the_threshold() {
        assert_eq!(spread_cell(20, 7, MIN_TEMPLATES_TO_SPREAD).len(), SPREAD * SPREAD);
    }

    #[test]
    fn spread_cell_leaves_sparse_cells_as_one_point() {
        assert_eq!(spread_cell(20, 7, MIN_TEMPLATES_TO_SPREAD - 1), vec![(20.0, 7.0, 3.0)]);
    }
}
