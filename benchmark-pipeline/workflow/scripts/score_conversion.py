#!/usr/bin/env python3
"""Score one tool's conversion-failure calls against holodeck's golden BAM.

Truth comes from the golden BAM's per-molecule ``cf:i`` tag (identical on every
record of a QNAME): ``cf=1`` means the molecule was a whole-molecule conversion
failure (the positive class), ``cf=0`` means it converted normally.

Each template lands in exactly one cell of a 2x3 contingency table:

                       called converted   called non-converted   called mixed
    truth converted    (cf=0)  TN                FP                  cf0_mixed
    truth non-conv.    (cf=1)  FN                TP                  cf1_mixed

"called" is read from the tool's own output, per a tool-specific rule:

  * methylsieve / neb : a record is flagged if it carries the verdict tag
                        XX:Z:UC or has the 0x200 QC-fail bit set.
  * biscuit / bismark : the tool *output* is itself the set of reads it called
                        unconverted (``biscuit bsconv -m N -v`` /
                        bismark ``*_removed_seqs.bam``), so every record present
                        is, by construction, a flagged mate.

A template is "called non-converted" when all its primary mates are flagged,
"converted" when none are, and "mixed" when only some are. methylsieve decides
per-template and propagates the call to every record, so its mixed count is
structurally 0; the per-read tools can produce real mixed templates, which the
table surfaces directly. The downstream mixed_policy (in collate.py) decides how
mixed folds into sensitivity/specificity.

Memory is bounded: only QNAMEs with at least one flagged mate are held in a
dict; the golden BAM is streamed in QNAME-adjacent groups (it is query-grouped).
"""

import argparse
import sys

import pysam


def is_primary(read):
    """A primary alignment: not secondary (0x100), not supplementary (0x800)."""
    return not (read.is_secondary or read.is_supplementary)


def _tagged_unconverted(read):
    """True if methylsieve/NEB flagged this record (XX:Z:UC tag or 0x200)."""
    if read.has_tag("XX") and str(read.get_tag("XX")) == "UC":
        return True
    return bool(read.flag & 0x200)


def flagged_mate_counts(tool_output, caller):
    """Return {qname: n_flagged_primary_records} for templates the tool flagged.

    Only QNAMEs with >= 1 flagged primary record are stored (bounded by the
    dataset's conversion-failure rate), so converted templates cost no memory.
    """
    if caller in ("methylsieve", "neb"):
        is_flagged = _tagged_unconverted
    elif caller in ("biscuit", "bismark"):
        # The tool's output IS the set of reads it called unconverted, so every
        # record present is flagged.
        def is_flagged(_read):
            return True
    else:
        sys.exit(f"unknown caller: {caller}")

    counts = {}
    with pysam.AlignmentFile(tool_output, check_sq=False) as fh:
        for read in fh:
            if not is_primary(read):
                continue
            if is_flagged(read):
                counts[read.query_name] = counts.get(read.query_name, 0) + 1
    return counts


def iter_templates(golden_bam):
    """Yield (qname, cf, n_primary_mates) per QNAME group of the golden BAM.

    The golden BAM is query-grouped (all records of a QNAME adjacent), so groups
    are runs of consecutive records sharing a name.
    """
    cur_name = None
    cf = None
    n_mates = 0
    with pysam.AlignmentFile(golden_bam, check_sq=False) as fh:
        for read in fh:
            if read.query_name != cur_name:
                if cur_name is not None:
                    yield cur_name, cf, n_mates
                cur_name = read.query_name
                cf = None
                n_mates = 0
            if cf is None and read.has_tag("cf"):
                cf = int(read.get_tag("cf"))
            if is_primary(read):
                n_mates += 1
    if cur_name is not None:
        yield cur_name, cf, n_mates


def score(golden_bam, tool_output, caller):
    flagged = flagged_mate_counts(tool_output, caller)

    # cells[cf][bucket]; bucket in {converted, nonconv, mixed}
    cells = {0: {"converted": 0, "nonconv": 0, "mixed": 0},
             1: {"converted": 0, "nonconv": 0, "mixed": 0}}
    n_no_truth = 0

    for qname, cf, n_mates in iter_templates(golden_bam):
        if cf is None or cf not in (0, 1):
            n_no_truth += 1
            continue
        n_flagged = flagged.get(qname, 0)
        if n_mates <= 0:
            bucket = "converted"  # no evidence records -> treated as converted
        elif n_flagged == 0:
            bucket = "converted"
        elif n_flagged >= n_mates:
            bucket = "nonconv"
        else:
            bucket = "mixed"
        cells[cf][bucket] += 1

    return cells, n_no_truth


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--golden", required=True, help="golden BAM with cf:i truth")
    ap.add_argument("--tool-output", required=True, help="the tool's output BAM/SAM")
    ap.add_argument("--caller", required=True,
                    choices=["methylsieve", "neb", "biscuit", "bismark"])
    ap.add_argument("--layout", required=True, choices=["PE", "SE"])
    ap.add_argument("--output", required=True, help="metrics TSV (single data row)")
    args = ap.parse_args()

    cells, n_no_truth = score(args.golden, args.tool_output, args.caller)

    n_templates = sum(sum(b.values()) for b in cells.values())
    row = {
        "caller": args.caller,
        "layout": args.layout,
        "cf0_called_converted": cells[0]["converted"],
        "cf0_called_nonconv":   cells[0]["nonconv"],
        "cf0_called_mixed":     cells[0]["mixed"],
        "cf1_called_converted": cells[1]["converted"],
        "cf1_called_nonconv":   cells[1]["nonconv"],
        "cf1_called_mixed":     cells[1]["mixed"],
        "n_truth_converted":    sum(cells[0].values()),
        "n_truth_nonconv":      sum(cells[1].values()),
        "n_templates":          n_templates,
        "n_no_truth":           n_no_truth,
    }
    cols = list(row.keys())
    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        out.write("\t".join(str(row[c]) for c in cols) + "\n")

    if n_no_truth:
        print(f"warning: {n_no_truth} templates lacked a cf tag and were skipped",
              file=sys.stderr)


if __name__ == "__main__":
    main()
