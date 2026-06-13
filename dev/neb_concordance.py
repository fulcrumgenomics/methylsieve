#!/usr/bin/env python3
"""Concordance harness: methylsieve vs NEB mark-nonconverted-reads.

Generates synthetic directional-EM-seq reads programmatically (no committed
fixtures), runs both tools on identical input, and compares the set of reads
tagged ``XX:Z:UC``.

Two divergences are EXPECTED and are characterized rather than treated as
failures:

  1. NEB requires ``is_proper_pair`` and so passes *single-end* reads through
     untouched; methylsieve evaluates them. (methylsieve is a superset here.)
  2. NEB decides per *read*; methylsieve decides per *template* (R1+R2 evidence
     summed). They agree only when each mate is independently decisive, so the
     concordance assertion runs on matched-class proper pairs; a deliberately
     discordant pair is reported separately.

Usage:
    neb-venv/bin/python dev/neb_concordance.py \
        --methylsieve target/release/methylsieve \
        --neb ../mark-nonconverted-reads/mark-nonconverted-reads.py
"""
import argparse
import os
import subprocess
import sys
import tempfile

import pysam

# Reference: a CpA region (top-strand C's, monitor C) followed by a
# CpA-equivalent region (bottom-strand G's, monitor G). Every monitored
# cytosine is in CpH context so methylsieve's default CpH == NEB's non-CpG.
REGION_TOP = "CA" * 20  # C at even offsets 0..38, each followed by A → CpA
REGION_BOT = "TG" * 20  # G at odd offsets 1..39, each preceded by T → CpA-equiv
REFERENCE = REGION_TOP + REGION_BOT
TOP_START = 0            # 0-based start of a 30 bp read in the top region
BOT_START = 40           # 0-based start of a 30 bp read in the bottom region
READ_LEN = 30
QUAL = "I" * READ_LEN    # Phred 40


def top_read(unconv):
    """A 30 bp top-strand read with `unconv` unconverted CpA cytosines."""
    bases = []
    for i in range(READ_LEN):
        if i % 2 == 1:
            bases.append("A")            # matches reference
        elif (i // 2) < unconv:
            bases.append("C")            # unconverted
        else:
            bases.append("T")            # converted
    return "".join(bases)


def bot_read(unconv):
    """A 30 bp bottom-strand read (ref-forward SEQ) with `unconv` unconverted Gs."""
    bases = []
    for i in range(READ_LEN):
        if i % 2 == 0:
            bases.append("T")            # matches reference
        elif (i // 2) < unconv:
            bases.append("G")            # unconverted
        else:
            bases.append("A")            # converted
    return "".join(bases)


# SAM flag building blocks.
PAIRED, PROPER, RUNMAP, MUNMAP = 0x1, 0x2, 0x4, 0x8
REV, MREV, R1, R2 = 0x10, 0x20, 0x40, 0x80


class Sam:
    def __init__(self):
        self.lines = [f"@SQ\tSN:chr1\tLN:{len(REFERENCE)}"]

    def rec(self, qname, flag, pos, seq):
        self.lines.append(
            f"{qname}\t{flag}\tchr1\t{pos}\t60\t{READ_LEN}M\t=\t{pos}\t0\t{seq}\t{QUAL}"
        )

    def ot_pair(self, qname, u_r1, u_r2):
        self.rec(qname, PAIRED | PROPER | R1 | MREV, TOP_START + 1, top_read(u_r1))
        self.rec(qname, PAIRED | PROPER | R2 | REV, TOP_START + 1, top_read(u_r2))

    def ob_pair(self, qname, u_r1, u_r2):
        self.rec(qname, PAIRED | PROPER | R1 | REV, BOT_START + 1, bot_read(u_r1))
        self.rec(qname, PAIRED | PROPER | R2 | MREV, BOT_START + 1, bot_read(u_r2))

    def se(self, qname, flag, seq, pos):
        self.rec(qname, flag, pos, seq)

    def write(self, path):
        with open(path, "w") as fh:
            fh.write("\n".join(self.lines) + "\n")


def read_key(read):
    if read.is_read1:
        end = "R1"
    elif read.is_read2:
        end = "R2"
    else:
        end = "SE"
    return f"{read.query_name}/{end}"


def tagged_set(path, mode):
    """Return the set of read keys carrying an XX tag in a SAM/BAM file."""
    out = set()
    with pysam.AlignmentFile(path, mode) as af:
        for read in af:
            if read.has_tag("XX"):
                out.add(read_key(read))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--methylsieve", default="target/release/methylsieve")
    ap.add_argument("--neb", default="../mark-nonconverted-reads/mark-nonconverted-reads.py")
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix="neb-concordance-")
    ref = os.path.join(work, "chr1.fa")
    with open(ref, "w") as fh:
        fh.write(f">chr1\n{REFERENCE}\n")
    subprocess.run(["samtools", "faidx", ref], check=True)

    sam = Sam()
    # Matched-class proper pairs (per-read == per-template) — concordance set.
    concordant = []
    for u in (0, 3, 4, 5, 15):
        sam.ot_pair(f"ot_u{u}", u, u)
        sam.ob_pair(f"ob_u{u}", u, u)
        for end in ("R1", "R2"):
            concordant.append(f"ot_u{u}/{end}")
            concordant.append(f"ob_u{u}/{end}")
    # Deliberately discordant pair (R1 unconverted, R2 fully converted).
    sam.ot_pair("ot_discordant", 5, 0)
    # Single-end reads (forward and reverse) with 5 unconverted each.
    sam.se("se_fwd", 0, top_read(5), TOP_START + 1)
    sam.se("se_rev", REV, bot_read(5), BOT_START + 1)

    in_sam = os.path.join(work, "in.sam")
    sam.write(in_sam)

    # Run methylsieve (BAM out). -q 0 to match NEB's lack of a base-quality
    # filter; default contexts (CpH), count 3, tag XX:Z:UC.
    ms_bam = os.path.join(work, "ms.bam")
    subprocess.run(
        [args.methylsieve, "-i", in_sam, "-o", ms_bam, "-r", ref,
         "--min-base-quality", "0", "-q"],
        check=True,
    )

    # Run NEB (SAM out).
    neb_sam = os.path.join(work, "neb.sam")
    subprocess.run(
        [sys.executable, args.neb, "--bam", in_sam, "--reference", ref, "--out", neb_sam],
        check=True,
    )

    ms = tagged_set(ms_bam, "rb")
    neb = tagged_set(neb_sam, "r")

    concordant = set(concordant)
    ms_c = ms & concordant
    neb_c = neb & concordant

    print(f"work dir: {work}")
    print(f"methylsieve tagged (all): {sorted(ms)}")
    print(f"NEB tagged (all):         {sorted(neb)}")
    print()
    print("=== Concordance on matched-class proper pairs ===")
    only_ms = ms_c - neb_c
    only_neb = neb_c - ms_c
    print(f"  methylsieve-only: {sorted(only_ms)}")
    print(f"  NEB-only:         {sorted(only_neb)}")
    ok = not only_ms and not only_neb

    print()
    print("=== Expected divergences (not failures) ===")
    se_ms = {k for k in ms if k.startswith("se_")}
    se_neb = {k for k in neb if k.startswith("se_")}
    print(f"  single-end tagged by methylsieve: {sorted(se_ms)}")
    print(f"  single-end tagged by NEB:         {sorted(se_neb)}  (NEB skips non-proper-pairs)")
    disc_ms = {k for k in ms if k.startswith("ot_discordant")}
    disc_neb = {k for k in neb if k.startswith("ot_discordant")}
    print(f"  discordant pair tagged by methylsieve: {sorted(disc_ms)}  (per-template)")
    print(f"  discordant pair tagged by NEB:         {sorted(disc_neb)}  (per-read)")

    # Validate the expected divergences too.
    se_ok = se_ms == {"se_fwd/SE", "se_rev/SE"} and se_neb == set()
    disc_ok = disc_ms == {"ot_discordant/R1", "ot_discordant/R2"} and disc_neb == {"ot_discordant/R1"}

    print()
    if ok and se_ok and disc_ok:
        print("CONCORDANT: methylsieve matches NEB on the shared code path; "
              "the two documented divergences behaved exactly as expected.")
        return 0
    print("MISMATCH — see above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
