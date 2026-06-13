#!/usr/bin/env python3
"""Join every per-(dataset x tool x rep) score + timing into one big TSV.

Reads each ``results/score/{id}/{tool}/rep{rep}/metrics.tsv`` (the 2x3
contingency from score_conversion.py) and its sibling
``results/run/{id}/{tool}/rep{rep}/time.txt`` (GNU ``time -v``), prepends the
dataset parameters and provenance (host, arch, tool versions), derives
sensitivity/specificity/PPV/F1 under the configured ``mixed_policy``, and writes
``results/benchmark.tsv`` — one row per (dataset, tool, rep).
"""

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_gnu_time import parse as parse_time  # noqa: E402

DATASET_COLS = ["layout", "read_length", "fragment_mean", "fragment_stddev",
                "coverage", "conversion_rate", "failure_rate",
                "min_error_rate", "max_error_rate", "seed"]
TIME_COLS = ["wall_s", "user_s", "sys_s", "cpu_percent", "max_rss_kb", "exit_status"]
CONTINGENCY_COLS = ["cf0_called_converted", "cf0_called_nonconv", "cf0_called_mixed",
                    "cf1_called_converted", "cf1_called_nonconv", "cf1_called_mixed"]


def load_datasets(path):
    with open(path) as fh:
        return {r["id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def read_metrics(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def fold_contingency(m, policy):
    """Fold the 6 cells into (TP, FN, FP, TN) under the mixed policy.

    Positive class = truth non-converted (cf=1). `mixed` templates (mates
    disagree) fold per policy: counted as non-converted (positive), as
    converted (negative), or excluded entirely (`separate`).
    """
    cf0c, cf0n, cf0m = (int(m[c]) for c in
                        ["cf0_called_converted", "cf0_called_nonconv", "cf0_called_mixed"])
    cf1c, cf1n, cf1m = (int(m[c]) for c in
                        ["cf1_called_converted", "cf1_called_nonconv", "cf1_called_mixed"])
    if policy == "non-converted":
        return cf1n + cf1m, cf1c, cf0n + cf0m, cf0c
    if policy == "converted":
        return cf1n, cf1c + cf1m, cf0n, cf0c + cf0m
    if policy == "separate":
        return cf1n, cf1c, cf0n, cf0c
    sys.exit(f"unknown mixed_policy: {policy!r} (use converted|non-converted|separate)")


def ratio(num, den):
    return round(num / den, 6) if den else ""


def per_second(count, seconds):
    """Throughput = count / wall seconds, rounded; '' if seconds missing or 0."""
    try:
        s = float(seconds)
    except (TypeError, ValueError):
        return ""
    return round(count / s) if s > 0 else ""


_version_cache = {}


def tool_version(tool, methylsieve_bin, holodeck_bin):
    if tool in _version_cache:
        return _version_cache[tool]
    cmds = {
        "methylsieve": [methylsieve_bin, "--version"],
        "biscuit": ["biscuit", "version"],
        "bismark": ["pixi", "run", "-e", "bismark", "bismark", "--version"],
        "neb": None,  # the NEB script exposes no version
    }
    cmd = cmds.get(tool)
    version = ""
    if cmd is not None:
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            for line in (out.stdout + "\n" + out.stderr).splitlines():
                if line.strip():
                    version = line.strip()
                    break
        except Exception:
            version = ""
    _version_cache[tool] = version
    return version


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--mixed-policy", required=True)
    ap.add_argument("--methylsieve-bin", required=True)
    ap.add_argument("--holodeck-bin", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    datasets = load_datasets(args.datasets)
    host = json.loads(Path(args.host).read_text())
    holodeck_ver = _capture([args.holodeck_bin, "--version"])
    results = Path(args.results_dir)

    rows = []
    for metrics_path in sorted(results.glob("score/*/*/rep*/metrics.tsv")):
        rep_dir = metrics_path.parent
        rep = int(rep_dir.name.removeprefix("rep"))
        tool = rep_dir.parent.name
        dataset_id = rep_dir.parent.parent.name

        m = read_metrics(metrics_path)
        if not m:
            print(f"warning: empty metrics {metrics_path}", file=sys.stderr)
            continue
        timing = parse_time(results / "run" / dataset_id / tool / f"rep{rep}" / "time.txt")
        ds = datasets.get(dataset_id, {})

        tp, fn, fp, tn = fold_contingency(m, args.mixed_policy)
        n_templates = int(m["n_templates"])
        n_mixed = int(m["cf0_called_mixed"]) + int(m["cf1_called_mixed"])

        row = {"id": dataset_id}
        row.update({c: ds.get(c, "") for c in DATASET_COLS})
        row["tool"] = tool
        row["rep"] = rep
        row.update({c: timing.get(c, "") for c in TIME_COLS})
        # Throughput from wall time. The golden BAM has exactly `mates` primary
        # records per template (2 for PE, 1 for SE), so reads = templates*mates.
        mates = 2 if (ds.get("layout") or m.get("layout") or "").upper() == "PE" else 1
        n_reads = n_templates * mates
        read_len = ds.get("read_length")
        wall = timing.get("wall_s", "")
        row["reads_per_s"] = per_second(n_reads, wall)
        row["templates_per_s"] = per_second(n_templates, wall)
        row["bases_per_s"] = per_second(n_reads * int(read_len), wall) if read_len else ""
        row.update({c: int(m[c]) for c in CONTINGENCY_COLS})
        row["n_templates"] = n_templates
        row["sensitivity"] = ratio(tp, tp + fn)
        row["specificity"] = ratio(tn, tn + fp)
        row["ppv"] = ratio(tp, tp + fp)
        row["f1"] = ratio(2 * tp, 2 * tp + fp + fn)
        row["mixed_rate"] = ratio(n_mixed, n_templates)
        row["mixed_policy"] = args.mixed_policy
        row["host"] = host.get("hostname", "")
        row["arch"] = host.get("arch", "")
        row["cpu_model"] = host.get("cpu_model", "")
        row["tool_version"] = tool_version(tool, args.methylsieve_bin, args.holodeck_bin)
        row["holodeck_version"] = holodeck_ver
        rows.append(row)

    if not rows:
        sys.exit("no metrics found under results/score/")

    df = pd.DataFrame(rows).sort_values(["id", "tool", "rep"]).reset_index(drop=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"wrote {len(df)} rows -> {args.output}")


def _capture(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        for line in (out.stdout + "\n" + out.stderr).splitlines():
            if line.strip():
                return line.strip()
    except Exception:
        pass
    return ""


if __name__ == "__main__":
    main()
