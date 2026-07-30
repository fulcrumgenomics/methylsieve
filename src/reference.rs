//! Reference FASTA loading and per-position context lookup.
//!
//! Ported from riker's `fasta.rs` / `sequence_dict.rs`. Uses `noodles-fasta`
//! only to open the file and parse the `.fai` index; the sequence bytes are
//! read in-house (see below).
//!
//! The reference is loaded **once at startup**, packed **2 bits per base**
//! (A/C/G/T → 0/1/2/3; every non-ACGT base — `N` and IUPAC ambiguity codes —
//! folds to A). The per-record hot path compares reference bases against read
//! bases — which [`fgumi_raw_bam::RawRecord::get_base`] returns as 4-bit BAM
//! codes (A=1, C=2, G=4, T=8) — by decoding the 2-bit value back to that code
//! space, or, on the hottest predicates, comparing directly in 2-bit space.
//!
//! Loading is **index-driven**: when a `.fai` is present each contig's byte
//! span is read in one bulk `read_exact` and newlines are stripped by the
//! known line geometry while packing — no per-line scanning, and a single pass
//! over the bases. Without a `.fai` we fall back to a sequential `noodles`
//! read; both paths produce identical packed output.

use std::collections::HashMap;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use anyhow::{Context as _, Result, anyhow, bail};
use noodles_fasta as fasta;
use noodles_sam::Header;

// ── 4-bit BAM base codes (subset we care about) ─────────────────────────────

/// 4-bit BAM code for cytosine.
pub(crate) const BASE_C: u8 = 2;
/// 4-bit BAM code for guanine.
pub(crate) const BASE_G: u8 = 4;
/// 4-bit BAM code for adenine.
pub(crate) const BASE_A: u8 = 1;
/// 4-bit BAM code for thymine.
pub(crate) const BASE_T: u8 = 8;
/// 4-bit BAM code for an unknown / ambiguous base. Completes the documented
/// code family; the 2-bit store never produces it (non-ACGT folds to A), so it
/// is only referenced by the test encoder.
#[allow(dead_code)]
pub(crate) const BASE_N: u8 = 15;

/// Sentinel in [`REF_ASCII_TO_2BIT`] for a byte that is not a valid sequence
/// character. The index-driven loader locates bases by `.fai` arithmetic, so a
/// corrupt index or malformed FASTA can slip a structural byte (`>`, a digit,
/// whitespace, a stray newline) into what we treat as sequence; folding those
/// silently to A would build a wrong reference, so the loader rejects them.
const INVALID_BASE: u8 = 0xFF;

/// ASCII reference base → 2-bit code (A=0, C=1, G=2, T=3). Any other **letter** —
/// IUPAC ambiguity codes and `N` — folds to A (0); this matches the prior
/// ASCII→4-bit→2-bit fold (non-ACGT became N then A): folding N→A never creates a
/// spurious monitored site (A is never the monitored C/G), and only gives a
/// former-N neighbor a concrete context instead of being skipped — immaterial in
/// practice (assembly gaps). Non-letter bytes map to [`INVALID_BASE`] so the
/// loader can reject them rather than silently fold structural bytes to A.
const REF_ASCII_TO_2BIT: [u8; 256] = build_ref_2bit();

const fn build_ref_2bit() -> [u8; 256] {
    let mut t = [INVALID_BASE; 256];
    // Every ASCII letter is a valid sequence character; non-ACGT (N / IUPAC)
    // folds to A (0). A/a are left at 0 by the ACGT assignments below.
    let mut b = b'A';
    while b <= b'Z' {
        t[b as usize] = 0;
        t[(b + 32) as usize] = 0; // lowercase
        b += 1;
    }
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}

// ── Conversion context ──────────────────────────────────────────────────────

/// The dinucleotide methylation context of a monitored cytosine, classified by
/// the base immediately 3' of the C (on the strand the molecule came from).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Context {
    /// C followed by A.
    CpA,
    /// C followed by C.
    CpC,
    /// C followed by G.
    CpG,
    /// C followed by T.
    CpT,
}

impl Context {
    /// All four contexts in stable index order (CpA, CpC, CpG, CpT).
    pub(crate) const ALL: [Context; 4] = [Context::CpA, Context::CpC, Context::CpG, Context::CpT];

    /// Stable index 0..4 for array-indexed per-context counters
    /// (CpA, CpC, CpG, CpT).
    #[inline]
    #[must_use]
    pub(crate) fn index(self) -> usize {
        match self {
            Context::CpA => 0,
            Context::CpC => 1,
            Context::CpG => 2,
            Context::CpT => 3,
        }
    }

    /// Canonical `CpX` label for TSV output.
    #[must_use]
    pub(crate) fn label(self) -> &'static str {
        match self {
            Context::CpA => "CpA",
            Context::CpC => "CpC",
            Context::CpG => "CpG",
            Context::CpT => "CpT",
        }
    }
}

// ── Reference access (2-bit packed) ──────────────────────────────────────────

/// Per-contig accessor over the 2-bit packed store, returning 4-bit BAM base
/// codes by position. Four bases per byte (~0.25 byte/base); only A/C/G/T are
/// representable, so all non-ACGT bases (N, IUPAC ambiguity) were folded to A at
/// load. This preserves monitored-site detection exactly (A is never the
/// monitored C/G) but means a monitored C/G adjacent to a former-N gets a
/// concrete context instead of being skipped (immaterial in practice — these
/// sit in assembly gaps).
///
/// The hot accessors ([`Self::for_each_monitored`], [`Self::ctx_top`],
/// [`Self::ctx_bottom`]) work directly in 2-bit space, so a scan never decodes to
/// the 4-bit code per base.
#[derive(Clone, Copy)]
pub(crate) struct TwoBitCodes<'a> {
    data: &'a [u8],
    len: usize,
}
impl TwoBitCodes<'_> {
    /// Number of bases in the contig. (No `is_empty`: real contigs are never
    /// empty and nothing branches on emptiness.)
    #[inline]
    #[allow(clippy::len_without_is_empty)]
    pub(crate) fn len(&self) -> usize {
        self.len
    }
    /// The 4-bit BAM base code at `pos` (`pos < len`). Production never decodes a
    /// whole base — it scans with [`Self::for_each_monitored`] and classifies with
    /// [`Self::ctx_top`] / [`Self::ctx_bottom`], all in 2-bit space. Tests use this.
    #[cfg(test)]
    #[inline]
    pub(crate) fn code(&self, pos: usize) -> u8 {
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        // 2-bit value → 4-bit BAM code without a table lookup: the codes for
        // A/C/G/T are 1/2/4/8 = 1 << (0/1/2/3).
        1u8 << val
    }
    /// Whether the base at `pos` equals the 4-bit BAM `code` (a monitored base,
    /// C=2 or G=4). Compares in 2-bit space: `code == 1 << val`, so
    /// `val == log2(code)`.
    ///
    /// The scan itself runs on [`Self::for_each_monitored`], which tests 32 bases
    /// at a time; this single-base form remains as the scalar definition that the
    /// bulk scan is asserted equal to.
    #[cfg(test)]
    #[inline]
    pub(crate) fn monitors(&self, pos: usize, code: u8) -> bool {
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        val == code.trailing_zeros() as u8
    }
    /// Call `f(pos)` for every position in `[start, end)` whose base equals the
    /// monitored 4-bit `code`, ascending. This is the tally's per-aligned-base
    /// reference scan.
    ///
    /// `end` is clamped to [`Self::len`], so a span past the contig simply finds
    /// nothing out there. Every caller already clips, so exceeding it is a caller
    /// bug and trips a debug assertion; the clamp is what keeps a release build
    /// from indexing past the packed store instead.
    ///
    /// Only ~21% of reference positions are the monitored C/G, so probing them
    /// one at a time spends most of its work rejecting. Instead this tests 32
    /// bases per 64-bit word: XOR against `code` splatted into every 2-bit field
    /// makes matching fields `00`, and a zero-field detect turns those into a
    /// bitmask iterated by `trailing_zeros`. The visited set is exactly
    /// `(start..end).filter(|p| monitors(p, code))` for the single-base predicate
    /// the tests assert this against.
    ///
    /// Both span edges are handled by clearing bits rather than by scalar loops,
    /// and the contig's final window is zero-filled where the packed store runs out
    /// (2-bit `0` is A, never a monitored base). That leaves **one** `f` call site,
    /// which is what lets the caller's per-site body inline into this loop — with
    /// three call sites LLVM outlines it and charges a call per site, undoing the
    /// win.
    #[inline]
    pub(crate) fn for_each_monitored(
        &self,
        start: usize,
        end: usize,
        code: u8,
        mut f: impl FnMut(usize),
    ) {
        // A span past the contig is a caller bug — every caller clips to the
        // contig — so say so loudly in debug. Clamp anyway: without it a release
        // build indexes past the packed store and panics obscurely instead of
        // simply finding no monitored bases out there.
        debug_assert!(end <= self.len, "span end {end} past contig length {}", self.len);
        let end = end.min(self.len);
        if start >= end {
            return;
        }
        // `code` as a 2-bit value repeated across all 32 fields. `0x5555… * v`
        // carries nothing for `v <= 3`, which holds for any single-base code.
        let splat = 0x5555_5555_5555_5555u64 * u64::from(code.trailing_zeros());

        // Windows are byte-aligned so a little-endian `u64` at byte `base / 4`
        // holds the base at `base + i` in bits `2i` — the match bitmask then
        // indexes positions directly. The first window may start below `start`.
        let mut base = start - start % 4;
        while base < end {
            let byte = base / 4;
            let word = match self.data.get(byte..byte + 8) {
                Some(chunk) => u64::from_le_bytes(chunk.try_into().expect("8 bytes")),
                // Final window of the contig: pad the missing bytes with zeros.
                None => {
                    let mut buf = [0u8; 8];
                    let avail = &self.data[byte..];
                    buf[..avail.len()].copy_from_slice(avail);
                    u64::from_le_bytes(buf)
                }
            };
            let zero = !(word ^ splat);
            let mut matches = zero & (zero >> 1) & 0x5555_5555_5555_5555;
            // Clear matches outside `[start, end)`. Each shift amount is < 64:
            // `start - base` is at most 3, and `end - base` is only used when it
            // is under 32.
            if base < start {
                matches &= !((1u64 << (2 * (start - base))) - 1);
            }
            if end - base < 32 {
                matches &= (1u64 << (2 * (end - base))) - 1;
            }
            while matches != 0 {
                f(base + matches.trailing_zeros() as usize / 2);
                matches &= matches - 1; // clear the lowest set bit
            }
            base += 32;
        }
    }
    /// Top-strand context implied by the base at `pos`, which the caller passes as
    /// the monitored C's 3' neighbor (`gp + 1`).
    #[inline]
    pub(crate) fn ctx_top(&self, pos: usize) -> Option<Context> {
        // 2-bit value 0/1/2/3 (A/C/G/T) indexes the contexts CpA/CpC/CpG/CpT
        // directly, skipping the `1 << val` decode and the `match`. There is no N
        // in 2-bit (folded to A), so the context is always defined.
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        Some(Context::ALL[val as usize])
    }
    /// Bottom-strand context implied by the base at `pos`, which the caller passes
    /// as the monitored G's 5' neighbor (`gp - 1`), complemented (A↔T, C↔G).
    #[inline]
    pub(crate) fn ctx_bottom(&self, pos: usize) -> Option<Context> {
        // Bottom strand uses the complement of the neighbor: A↔T, C↔G, which on
        // the 2-bit value is `3 - val` (A0↔T3, C1↔G2).
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        Some(Context::ALL[3 - val as usize])
    }
}

/// A packed contig: the 2-bit packed bytes plus the base count (which the
/// packing can't recover from the byte length alone). Held inside [`Reference`];
/// fields are private (access is via [`TwoBitCodes`]).
pub(crate) struct PackedContig {
    data: Vec<u8>,
    len: usize,
}

/// Pack a contiguous run of ASCII reference bases into the 2-bit store, skipping
/// line terminators by the FASTA's known geometry — no newline scanning.
///
/// `raw` is the on-disk byte span of one contig (sequence bytes plus line
/// terminators). `seq_len` is the base count; `line_bases` the bases per line
/// and `line_width` the bytes per line including the terminator (both from the
/// `.fai`). Bases are read `line_bases` at a time and the `line_width -
/// line_bases` terminator is skipped between lines. ASCII→2-bit translation
/// ([`REF_ASCII_TO_2BIT`]) is fused into this single pass, so the bases are
/// never materialized unpacked.
///
/// # Errors
/// Returns an error on a non-letter byte in the sequence span (see
/// [`INVALID_BASE`]) — a signal that the `.fai` geometry or FASTA is corrupt and
/// we would otherwise silently pack a wrong reference.
fn pack_twobit_from_lines(
    raw: &[u8],
    seq_len: usize,
    line_bases: usize,
    line_width: usize,
) -> Result<PackedContig> {
    let mut data = vec![0u8; seq_len.div_ceil(4)];
    if seq_len == 0 || line_bases == 0 {
        return Ok(PackedContig { data, len: seq_len });
    }
    let terminator_len = line_width.saturating_sub(line_bases);
    // Build each output byte from its four bases in a register (`acc`) and store
    // it once, rather than read-modify-writing `data` per base. Four consecutive
    // bases share one output byte, so the per-base RMW was four loads+ORs+stores
    // to the same byte on a store→load dependency chain; accumulating in a
    // register drops that chain and 3/4 of the bounds checks on `data`. The
    // partial byte (`acc`/`filled_bits`) carries across line boundaries, since a
    // line's base count need not be a multiple of four. 2-bit values are
    // little-endian within the byte (base 4k in bits 0-1), matching the unpacking
    // in `TwoBitCodes`.
    let mut acc: u8 = 0;
    let mut filled_bits: u8 = 0; // 0, 2, 4, or 6
    let mut out = 0; // next output byte index
    let mut src = 0; // byte index into `raw`
    let mut i = 0; // base index into the logical sequence
    while i < seq_len {
        let n = (seq_len - i).min(line_bases).min(raw.len().saturating_sub(src));
        if n == 0 {
            break; // malformed geometry / truncated span — stop rather than spin
        }
        // Slice once per line so the inner loop has no per-base bounds check on `raw`.
        for &b in &raw[src..src + n] {
            let code = REF_ASCII_TO_2BIT[b as usize];
            if code == INVALID_BASE {
                bail!(
                    "non-sequence byte {:#04x} ({:?}) where a reference base was expected — the \
                     .fai index or FASTA is likely corrupt",
                    b,
                    b as char
                );
            }
            acc |= code << filled_bits;
            filled_bits += 2;
            if filled_bits == 8 {
                data[out] = acc;
                out += 1;
                acc = 0;
                filled_bits = 0;
            }
        }
        i += n;
        src += n;
        if i < seq_len {
            src += terminator_len; // skip this line's terminator
        }
    }
    if filled_bits > 0 {
        data[out] = acc; // flush the final partial byte
    }
    Ok(PackedContig { data, len: seq_len })
}

/// The total on-disk byte span of a contig's sequence: complete lines (each
/// `line_width` bytes) plus the final partial line, which has no trailing
/// terminator. Robust to a missing final newline at EOF.
fn total_sequence_bytes(seq_len: usize, line_bases: usize, line_width: usize) -> usize {
    if seq_len == 0 || line_bases == 0 {
        return 0;
    }
    let complete = seq_len / line_bases;
    let remainder = seq_len % line_bases;
    if remainder > 0 {
        complete * line_width + remainder
    } else {
        // Exactly fills `complete` lines; the last has no trailing terminator.
        (complete - 1) * line_width + line_bases
    }
}

// ── Reference ───────────────────────────────────────────────────────────────

/// One contig's `.fai` geometry, copied out of the index so the index borrow
/// can be released before we seek/read the file (the bulk read borrows the
/// reader mutably).
struct FaiGeom {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

/// The full reference: every input `@SQ` contig, 2-bit packed, indexed by BAM
/// `tid` (the 0-based position of the `@SQ` line).
pub(crate) struct Reference {
    contigs: Vec<PackedContig>,
}

impl Reference {
    /// Load every contig named by `header`'s `@SQ` lines from the FASTA at
    /// `path`, in BAM tid order, 2-bit packed.
    ///
    /// When a `.fai` index sits beside the FASTA each contig is read by its
    /// byte span in one bulk `read_exact` ([`Self::load_indexed`]); without one
    /// we fall back to a sequential read ([`Self::load_sequential`]). Both paths
    /// produce identical packed output.
    ///
    /// Semantics:
    /// * The FASTA may be a **superset** of the BAM's contigs — extra FASTA
    ///   contigs are ignored.
    /// * A `@SQ` with no corresponding FASTA entry is a fatal error.
    /// * A `@SQ` whose length disagrees with the FASTA is a fatal error (the BAM
    ///   was aligned against a different reference).
    ///
    /// # Errors
    /// Returns an error if the FASTA / `.fai` cannot be opened, a contig is
    /// missing, a length mismatches, or a contig cannot be read.
    pub(crate) fn load(path: &Path, header: &Header) -> Result<Self> {
        if fai_path(path).is_some() {
            Self::load_indexed(path, header)
        } else {
            log::info!(
                "No .fai index beside {}; falling back to a sequential FASTA read \
                 (run `samtools faidx` for faster startup).",
                path.display()
            );
            Self::load_sequential(path, header)
        }
    }

    /// Index-driven load: bulk-read each contig's byte span and strip newlines
    /// by the `.fai` line geometry while packing. Peak load memory is one
    /// contig's raw bytes plus the (4×-smaller) accumulating packed store.
    fn load_indexed(path: &Path, header: &Header) -> Result<Self> {
        let mut reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(path).with_context(
                || format!("opening indexed FASTA {} (is there a .fai?)", path.display()),
            )?;

        // Copy the geometry out of the index so its borrow ends before the
        // mutable seek/read below.
        let geom: HashMap<String, FaiGeom> = reader
            .index()
            .as_ref()
            .iter()
            .map(|rec| {
                (
                    String::from_utf8_lossy(rec.name().as_ref()).into_owned(),
                    FaiGeom {
                        length: rec.length(),
                        offset: rec.offset(),
                        line_bases: rec.line_bases(),
                        line_width: rec.line_width(),
                    },
                )
            })
            .collect();

        let mut contigs = Vec::with_capacity(header.reference_sequences().len());
        for (name, map) in header.reference_sequences() {
            let name_str = std::str::from_utf8(name.as_ref())
                .map_err(|_| anyhow!("BAM @SQ name is not valid UTF-8"))?;
            let bam_len = usize::from(map.length());

            let g = geom.get(name_str).ok_or_else(|| {
                anyhow!(
                    "BAM contig '{name_str}' is not present in the reference FASTA index ({}). \
                     Every @SQ contig must exist in the reference.",
                    path.display()
                )
            })?;
            if g.length as usize != bam_len {
                bail!(
                    "Length mismatch for contig '{name_str}': BAM @SQ says {bam_len} bp but the \
                     reference FASTA says {} bp. The BAM was aligned against a different reference.",
                    g.length
                );
            }

            let (line_bases, line_width) = (g.line_bases as usize, g.line_width as usize);
            // A well-formed .fai has line_width >= line_bases > 0 (the terminator
            // is the difference). A malformed/hand-edited index that violates this
            // would make total_sequence_bytes mis-stride and silently pack a wrong
            // reference — reject it. noodles' .fai reader does not validate this.
            if bam_len > 0 && (line_bases == 0 || line_width < line_bases) {
                bail!(
                    "Corrupt .fai for contig '{name_str}': line_bases={line_bases}, \
                     line_width={line_width} (need line_width >= line_bases > 0)."
                );
            }
            let total = total_sequence_bytes(bam_len, line_bases, line_width);
            reader
                .get_mut()
                .seek(SeekFrom::Start(g.offset))
                .with_context(|| format!("seeking to contig '{name_str}' in {}", path.display()))?;
            let mut raw = vec![0u8; total];
            reader
                .get_mut()
                .read_exact(&mut raw)
                .with_context(|| format!("reading contig '{name_str}' from {}", path.display()))?;
            contigs.push(
                pack_twobit_from_lines(&raw, bam_len, line_bases, line_width)
                    .with_context(|| format!("packing contig '{name_str}'"))?,
            );
        }
        Ok(Reference { contigs })
    }

    /// Sequential fallback for an unindexed FASTA: read every record (noodles
    /// hands back newline-stripped bases), pack it, then assemble in `@SQ`
    /// order. The whole unpacked genome is not held — each record is packed and
    /// dropped as it is read.
    fn load_sequential(path: &Path, header: &Header) -> Result<Self> {
        let mut reader = fasta::io::reader::Builder
            .build_from_path(path)
            .with_context(|| format!("opening FASTA {}", path.display()))?;

        let mut by_name: HashMap<String, PackedContig> = HashMap::new();
        for result in reader.records() {
            let record = result
                .with_context(|| format!("reading a FASTA record from {}", path.display()))?;
            let name = String::from_utf8_lossy(record.name()).into_owned();
            let bases = record.sequence().as_ref();
            // Treat the stripped sequence as a single line (no terminators).
            let packed = pack_twobit_from_lines(bases, bases.len(), bases.len(), bases.len())
                .with_context(|| format!("packing contig '{name}'"))?;
            by_name.insert(name, packed);
        }

        let mut contigs = Vec::with_capacity(header.reference_sequences().len());
        for (name, map) in header.reference_sequences() {
            let name_str = std::str::from_utf8(name.as_ref())
                .map_err(|_| anyhow!("BAM @SQ name is not valid UTF-8"))?;
            let bam_len = usize::from(map.length());

            let packed = by_name.remove(name_str).ok_or_else(|| {
                anyhow!(
                    "BAM contig '{name_str}' is not present in the reference FASTA ({}). \
                     Every @SQ contig must exist in the reference.",
                    path.display()
                )
            })?;
            if packed.len != bam_len {
                bail!(
                    "Length mismatch for contig '{name_str}': BAM @SQ says {bam_len} bp but the \
                     reference FASTA says {} bp. The BAM was aligned against a different reference.",
                    packed.len
                );
            }
            contigs.push(packed);
        }
        Ok(Reference { contigs })
    }

    /// 2-bit-packed contig accessor for `tid`. `None` for an unmapped record
    /// (`tid < 0`) or a `tid` past the loaded contigs.
    #[inline]
    #[must_use]
    pub(crate) fn codes(&self, tid: i32) -> Option<TwoBitCodes<'_>> {
        if tid < 0 {
            return None;
        }
        self.contigs.get(tid as usize).map(|c| TwoBitCodes { data: &c.data, len: c.len })
    }

    /// Construct directly from already-encoded (4-bit BAM code) contigs, packed
    /// to the 2-bit store (test helper).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn from_encoded_contigs(contigs: Vec<Vec<u8>>) -> Self {
        Reference { contigs: contigs.iter().map(|c| pack_codes_twobit(c)).collect() }
    }
}

/// The sibling `.fai` path (`<fasta>.fai`) if it exists. Built by appending to
/// the path's raw bytes (not `Display`) so non-UTF-8 paths are preserved.
fn fai_path(path: &Path) -> Option<std::path::PathBuf> {
    let mut candidate = path.as_os_str().to_os_string();
    candidate.push(".fai");
    let candidate = std::path::PathBuf::from(candidate);
    candidate.exists().then_some(candidate)
}

/// Encode an ASCII reference base to its 4-bit BAM code (test-only convenience
/// for building encoded contigs and expected values; production loads pack
/// straight from ASCII to 2-bit via [`REF_ASCII_TO_2BIT`]).
#[cfg(test)]
#[must_use]
pub(crate) fn encode_ref_base(ascii: u8) -> u8 {
    match ascii {
        b'A' | b'a' => BASE_A,
        b'C' | b'c' => BASE_C,
        b'G' | b'g' => BASE_G,
        b'T' | b't' => BASE_T,
        _ => BASE_N,
    }
}

/// Pack 4-bit BAM codes into the 2-bit store (test-only; production packs
/// straight from ASCII via [`pack_twobit_from_lines`]). Non-ACGT folds to A.
#[cfg(test)]
fn pack_codes_twobit(codes: &[u8]) -> PackedContig {
    let mut data = vec![0u8; codes.len().div_ceil(4)];
    for (i, &c) in codes.iter().enumerate() {
        let val = match c {
            BASE_C => 1,
            BASE_G => 2,
            BASE_T => 3,
            _ => 0, // A and any non-ACGT
        };
        data[i >> 2] |= val << ((i & 3) * 2);
    }
    PackedContig { data, len: codes.len() }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;
    use crate::sam_reader::SamReader;

    #[test]
    fn ascii_to_4bit_codes() {
        assert_eq!(encode_ref_base(b'C'), BASE_C);
        assert_eq!(encode_ref_base(b'c'), BASE_C);
        assert_eq!(encode_ref_base(b'G'), BASE_G);
        assert_eq!(encode_ref_base(b'A'), BASE_A);
        assert_eq!(encode_ref_base(b'T'), BASE_T);
        assert_eq!(encode_ref_base(b'N'), BASE_N);
        assert_eq!(encode_ref_base(b'R'), BASE_N); // ambiguity → N
    }

    #[test]
    fn twobit_ctx_top_buckets_by_next_base() {
        // ref[i+1] = A/C/G/T ⇒ CpA/CpC/CpG/CpT.
        let codes: Vec<u8> = "ACGT".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        assert_eq!(view.ctx_top(0), Some(Context::CpA));
        assert_eq!(view.ctx_top(1), Some(Context::CpC));
        assert_eq!(view.ctx_top(2), Some(Context::CpG));
        assert_eq!(view.ctx_top(3), Some(Context::CpT));
    }

    #[test]
    fn twobit_ctx_bottom_uses_complement_of_neighbor() {
        // Bottom strand takes the complement of ref[i-1]: A→T, C→G, G→C, T→A,
        // so contexts come out CpT/CpG/CpC/CpA for neighbor A/C/G/T.
        let codes: Vec<u8> = "ACGT".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        assert_eq!(view.ctx_bottom(0), Some(Context::CpT)); // A
        assert_eq!(view.ctx_bottom(1), Some(Context::CpG)); // C
        assert_eq!(view.ctx_bottom(2), Some(Context::CpC)); // G
        assert_eq!(view.ctx_bottom(3), Some(Context::CpA)); // T
    }

    #[test]
    fn twobit_packing_round_trips_acgt_and_folds_n_to_a() {
        let codes: Vec<u8> = "CAGTNCCGTA".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        assert_eq!(view.len(), codes.len());
        for (i, &c) in codes.iter().enumerate() {
            let expected = if c == BASE_N { BASE_A } else { c }; // N folds to A
            assert_eq!(view.code(i), expected, "2-bit mismatch at {i}");
        }
        assert_eq!(packed.data.len(), codes.len().div_ceil(4));
        // The folded N must never read as the monitored C/G.
        assert_ne!(view.code(4), BASE_C);
        assert_ne!(view.code(4), BASE_G);
    }

    #[test]
    fn twobit_monitors_matches_code_compare() {
        // `monitors` (native 2-bit compare) must agree with `code() == code` for
        // every position and both monitored bases.
        let codes: Vec<u8> = "CACGCATTGCGNCAGTACG".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        for i in 0..view.len() {
            for &code in &[BASE_C, BASE_G] {
                assert_eq!(view.monitors(i, code), view.code(i) == code, "pos {i} code {code}");
            }
        }
    }

    /// Positions `for_each_monitored` visits, as a vector.
    fn monitored_positions(
        view: &TwoBitCodes<'_>,
        start: usize,
        end: usize,
        code: u8,
    ) -> Vec<usize> {
        let mut got = Vec::new();
        view.for_each_monitored(start, end, code, |p| got.push(p));
        got
    }

    /// The scalar predicate applied over the same span — the definition
    /// `for_each_monitored` must reproduce exactly.
    fn scalar_positions(view: &TwoBitCodes<'_>, start: usize, end: usize, code: u8) -> Vec<usize> {
        (start..end).filter(|&p| view.monitors(p, code)).collect()
    }

    #[test]
    fn for_each_monitored_matches_scalar_probe_at_every_offset_and_length() {
        // 200 bases with C/G at irregular positions, so matches land on both
        // parities and straddle byte (4-base) and word (32-base) boundaries.
        // Every start offset covers all four byte alignments of the head loop,
        // and lengths run out to the contig end so the partial-final-word mask
        // and the zero-filled final window are both exercised.
        let codes: Vec<u8> =
            (0..200u32).map(|i| encode_ref_base(b"ACGTGCTAAC"[(i as usize * 7) % 10])).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };

        for &code in &[BASE_C, BASE_G] {
            for start in 0..40 {
                for end in start..=view.len() {
                    assert_eq!(
                        monitored_positions(&view, start, end, code),
                        scalar_positions(&view, start, end, code),
                        "code={code} start={start} end={end}"
                    );
                }
            }
        }
    }

    #[test]
    fn for_each_monitored_handles_contigs_shorter_than_one_word() {
        // Under 32 bases no whole 64-bit word is loadable, so every span is a single
        // partial window whose tail is zero-filled past the packed store.
        for len in 1..40usize {
            let codes: Vec<u8> = (0..len).map(|i| encode_ref_base(b"CGAT"[(i * 3) % 4])).collect();
            let packed = pack_codes_twobit(&codes);
            let view = TwoBitCodes { data: &packed.data, len: packed.len };
            for &code in &[BASE_C, BASE_G] {
                assert_eq!(
                    monitored_positions(&view, 0, len, code),
                    scalar_positions(&view, 0, len, code),
                    "len={len} code={code}"
                );
            }
        }
    }

    #[test]
    #[should_panic(expected = "past contig length")]
    fn for_each_monitored_rejects_a_span_past_the_contig() {
        // Callers clip their span to the contig, so overrunning it is a caller bug
        // and debug builds say so. Release builds clamp instead, which is what keeps
        // the zero-fill window from indexing past the packed store.
        let codes: Vec<u8> = "CACACACACA".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        view.for_each_monitored(0, packed.len + 64, BASE_C, |_| {});
    }

    #[test]
    fn for_each_monitored_visits_nothing_for_an_empty_span() {
        let codes: Vec<u8> = "CCCCCCCCCC".bytes().map(encode_ref_base).collect();
        let packed = pack_codes_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        assert!(monitored_positions(&view, 5, 5, BASE_C).is_empty());
        // An all-C contig: every position in the span matches, none outside it.
        assert_eq!(monitored_positions(&view, 3, 7, BASE_C), vec![3, 4, 5, 6]);
        assert!(monitored_positions(&view, 0, 10, BASE_G).is_empty());
    }

    #[test]
    fn context_index_is_stable() {
        assert_eq!(Context::CpA.index(), 0);
        assert_eq!(Context::CpC.index(), 1);
        assert_eq!(Context::CpG.index(), 2);
        assert_eq!(Context::CpT.index(), 3);
    }

    // ── Fused index-driven loader ─────────────────────────────────────────────

    /// Build a `Header` with the given `(name, length)` `@SQ` lines.
    fn header_for(contigs: &[(&str, usize)]) -> Header {
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n");
        for (name, len) in contigs {
            sam.push_str(&format!("@SQ\tSN:{name}\tLN:{len}\n"));
        }
        let boxed: Box<dyn std::io::BufRead> = Box::new(std::io::Cursor::new(sam.into_bytes()));
        SamReader::new(boxed).read_header().unwrap()
    }

    /// Write a FASTA wrapping each sequence at `width` bp/line, computing and
    /// (optionally) writing the matching `.fai` from the byte layout we control.
    /// Returns the FASTA path.
    fn write_fasta(
        dir: &std::path::Path,
        seqs: &[(&str, &str)],
        width: usize,
        write_fai: bool,
    ) -> std::path::PathBuf {
        let fa = dir.join("ref.fasta");
        let mut content: Vec<u8> = Vec::new();
        let mut fai = String::new();
        for (name, seq) in seqs {
            content.extend_from_slice(format!(">{name}\n").as_bytes());
            let offset = content.len();
            let bytes = seq.as_bytes();
            let mut i = 0;
            while i < bytes.len() {
                let end = (i + width).min(bytes.len());
                content.extend_from_slice(&bytes[i..end]);
                content.push(b'\n');
                i = end;
            }
            // fai columns: NAME LENGTH OFFSET LINEBASES LINEWIDTH.
            fai.push_str(&format!("{name}\t{}\t{offset}\t{width}\t{}\n", seq.len(), width + 1));
        }
        std::fs::File::create(&fa).unwrap().write_all(&content).unwrap();
        if write_fai {
            std::fs::File::create(format!("{}.fai", fa.display()))
                .unwrap()
                .write_all(fai.as_bytes())
                .unwrap();
        }
        fa
    }

    /// All positions of every contig, as decoded 4-bit codes — the observable
    /// behavior, regardless of internal packing.
    fn decode_all(reference: &Reference, contigs: &[(&str, usize)]) -> Vec<Vec<u8>> {
        (0..contigs.len() as i32)
            .map(|tid| {
                let c = reference.codes(tid).unwrap();
                (0..c.len()).map(|p| c.code(p)).collect()
            })
            .collect()
    }

    /// A mixed payload: a multi-line contig with lowercase, `N`, and IUPAC
    /// ambiguity, length not a multiple of the line width (exercises the partial
    /// last line); plus a single-line contig (the `seq_len <= line_bases` path).
    const SEQS: &[(&str, &str)] = &[("chr1", "ACGTacgtNNRYACGTACGTAC"), ("chr2", "ACGT")];

    #[test]
    fn indexed_and_unindexed_load_produce_identical_codes() {
        let contigs: Vec<(&str, usize)> = SEQS.iter().map(|(n, s)| (*n, s.len())).collect();
        let header = header_for(&contigs);

        let idx_dir = tempfile::tempdir().unwrap();
        let idx_fa = write_fasta(idx_dir.path(), SEQS, 10, true);
        let indexed = Reference::load(&idx_fa, &header).unwrap();
        assert!(fai_path(&idx_fa).is_some(), "indexed case must have a .fai");

        let seq_dir = tempfile::tempdir().unwrap();
        let seq_fa = write_fasta(seq_dir.path(), SEQS, 10, false);
        let sequential = Reference::load(&seq_fa, &header).unwrap();
        assert!(fai_path(&seq_fa).is_none(), "unindexed case must have no .fai");

        assert_eq!(decode_all(&indexed, &contigs), decode_all(&sequential, &contigs));
    }

    #[test]
    fn loaded_codes_match_expected_bases_and_fold_non_acgt() {
        let contigs: Vec<(&str, usize)> = SEQS.iter().map(|(n, s)| (*n, s.len())).collect();
        let header = header_for(&contigs);
        let dir = tempfile::tempdir().unwrap();
        let fa = write_fasta(dir.path(), SEQS, 10, true);
        let reference = Reference::load(&fa, &header).unwrap();

        // chr1 = ACGTacgtNNRYACGTACGTAC: case-insensitive A/C/G/T, and N/R/Y → A.
        let expected: Vec<u8> = "ACGTACGTAAAAACGTACGTAC".bytes().map(encode_ref_base).collect();
        let chr1 = reference.codes(0).unwrap();
        assert_eq!(chr1.len(), expected.len());
        for (p, &e) in expected.iter().enumerate() {
            assert_eq!(chr1.code(p), e, "chr1 pos {p}");
        }
        // chr2 single-line contig reads back exactly.
        let chr2 = reference.codes(1).unwrap();
        assert_eq!(
            (0..chr2.len()).map(|p| chr2.code(p)).collect::<Vec<_>>(),
            vec![BASE_A, BASE_C, BASE_G, BASE_T]
        );
    }

    #[test]
    fn load_rejects_contig_length_mismatch() {
        // BAM @SQ says chr1 is longer than the FASTA actually provides.
        let header = header_for(&[("chr1", 999)]);
        let dir = tempfile::tempdir().unwrap();
        let fa = write_fasta(dir.path(), &[("chr1", "ACGTACGT")], 10, true);
        let err = Reference::load(&fa, &header).err().expect("should error").to_string();
        assert!(err.contains("Length mismatch"), "unexpected error: {err}");
    }

    #[test]
    fn load_rejects_missing_contig() {
        let header = header_for(&[("chrX", 4)]);
        let dir = tempfile::tempdir().unwrap();
        let fa = write_fasta(dir.path(), &[("chr1", "ACGT")], 10, true);
        let err = Reference::load(&fa, &header).err().expect("should error").to_string();
        assert!(err.contains("not present in the reference"), "unexpected error: {err}");
    }

    #[test]
    fn load_rejects_non_sequence_bytes() {
        // A structural/garbage byte (here a digit) in the sequence must be an
        // error, not silently folded to A — the arithmetic loader would otherwise
        // build a wrong reference from a corrupt .fai or malformed FASTA. IUPAC
        // letters (N/R/Y) are fine (see the fold tests); only non-letters bail.
        let header = header_for(&[("chr1", 10)]);
        let dir = tempfile::tempdir().unwrap();
        let fa = write_fasta(dir.path(), &[("chr1", "ACGT1CGTAC")], 10, true);
        let err = format!("{:#}", Reference::load(&fa, &header).err().expect("should error"));
        assert!(err.contains("non-sequence byte"), "unexpected error: {err}");
    }

    #[test]
    fn total_sequence_bytes_handles_partial_and_full_last_line() {
        // 22 bp at 10 bp/line (width 11 incl. \n): two full lines + "AC".
        assert_eq!(total_sequence_bytes(22, 10, 11), 11 + 11 + 2);
        // Exactly two full lines (20 bp): the last line has no trailing \n.
        assert_eq!(total_sequence_bytes(20, 10, 11), 11 + 10);
        // Single short line.
        assert_eq!(total_sequence_bytes(4, 10, 11), 4);
        assert_eq!(total_sequence_bytes(0, 10, 11), 0);
    }
}
