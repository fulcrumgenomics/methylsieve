#!/usr/bin/env python3
"""Render the README throughput chart from ``results/benchmark.tsv``.

Reads the collated benchmark TSV, averages the per-replicate throughput
(``reads_per_s``) for each ``(dataset, tool)`` pair, and draws a grouped bar
chart — one group per dataset, one bar per tool, with methylsieve highlighted.

Usage:
    python3 workflow/scripts/plot_benchmarks.py \\
        results/benchmark.tsv ../docs/img/throughput.png
"""
import csv
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Datasets in display order (paired-end first, then single-end), each with a
# compact two-line label for the x-axis.
DATASETS = [
    ("ins350_2x150", "350 bp insert\nPE  2×150"),
    ("cfdna_2x150", "cfDNA\nPE  2×150"),
    ("ins350_2x150_hifail", "350 bp insert\nPE  2×150 · 20% fail"),
    ("ins350_1x150", "350 bp insert\nSE  1×150"),
    ("cfdna-2nuc_1x400", "cfDNA di-nucleosome\nSE  1×400"),
]

# Tool order and palette: methylsieve gets the saturated brand colour; the
# peer tools use a muted, professional set so the eye lands on methylsieve.
TOOLS = [
    ("methylsieve", "methylsieve", "#0F9E8E"),
    ("biscuit", "biscuit", "#5B7FB4"),
    ("neb", "NEB", "#D69A3F"),
    ("bismark", "bismark", "#9AA0A6"),
]


def load_throughput(tsv_path):
    """Return ``{(dataset, tool): mean reads/s}`` averaged over replicates."""
    samples = defaultdict(list)
    with open(tsv_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            samples[(row["id"], row["tool"])].append(float(row["reads_per_s"]))
    return {key: float(np.mean(vals)) for key, vals in samples.items()}


def main():
    tsv_path = sys.argv[1] if len(sys.argv) > 1 else "results/benchmark.tsv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "throughput.png"
    throughput = load_throughput(tsv_path)

    x = np.arange(len(DATASETS))
    bar_w = 0.80 / len(TOOLS)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(figsize=(11.5, 5.6), dpi=200)

    for i, (tool, label, colour) in enumerate(TOOLS):
        heights = [throughput.get((ds, tool), 0.0) / 1000.0 for ds, _ in DATASETS]
        offsets = x + (i - (len(TOOLS) - 1) / 2) * bar_w
        bars = ax.bar(
            offsets, heights, bar_w, label=label,
            color=colour, edgecolor="white", linewidth=0.6, zorder=3,
        )
        for bar, height in zip(bars, heights):
            ax.annotate(
                f"{height:.0f}",
                (bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 2), textcoords="offset points",
                ha="center", va="bottom", fontsize=7.5, color="#333333",
            )

    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in DATASETS], fontsize=9.5)
    ax.set_ylabel("Throughput  (thousand reads / s)", fontsize=11.5)
    ax.set_title(
        "Conversion-failure calling throughput — higher is better",
        fontsize=13.5, fontweight="bold", pad=34,
    )
    ax.margins(x=0.02)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#DDDDDD", linewidth=0.8)
    ax.tick_params(axis="both", length=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#BBBBBB")

    ax.legend(
        ncol=4, frameon=False, loc="lower center",
        bbox_to_anchor=(0.5, 1.0), fontsize=10.5, columnspacing=1.8,
        handlelength=1.3,
    )
    fig.text(
        0.5, 0.01,
        "reads/s, mean of 2 replicates · 32× coverage of human chr21 · Apple M2 "
        "· biscuit input pre-sorted (sort untimed)",
        ha="center", fontsize=7.5, color="#888888",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
