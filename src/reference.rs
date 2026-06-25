//! Reference FASTA loading and per-position context lookup.
//!
//! Ported from riker's `fasta.rs` / `sequence_dict.rs` and adapted from the
//! umbrella `noodles` crate to the standalone `noodles-fasta` / `noodles-core`
//! crates so it shares one `noodles-core` with `noodles-sam` (see Cargo.toml).
//!
//! The reference is loaded **once at startup**. Every contig named by an input
//! `@SQ` line is read into a contiguous byte array, **encoded to the 4-bit BAM
//! base codes** (A=1, C=2, G=4, T=8, everything else=N=15) so the per-record
//! hot path compares reference bases against read bases — which
//! [`fgumi_raw_bam::RawRecord::get_base`] also returns as 4-bit codes — with a
//! plain nibble equality.

use std::collections::HashMap;
use std::io::{Seek, SeekFrom};

use anyhow::{Context as _, Result, anyhow, bail};
use noodles_core::Region;
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
/// 4-bit BAM code for an unknown / ambiguous base.
pub(crate) const BASE_N: u8 = 15;

/// ASCII reference base → 4-bit BAM code. Anything that isn't a plain A/C/G/T
/// (including IUPAC ambiguity codes and `N`) maps to `N` (15); such positions
/// are never a monitored C/G and never form a usable context, so they drop out
/// of tallying naturally.
const REF_ASCII_TO_CODE: [u8; 256] = build_ref_codes();

const fn build_ref_codes() -> [u8; 256] {
    let mut t = [BASE_N; 256];
    t[b'A' as usize] = BASE_A;
    t[b'a' as usize] = BASE_A;
    t[b'C' as usize] = BASE_C;
    t[b'c' as usize] = BASE_C;
    t[b'G' as usize] = BASE_G;
    t[b'g' as usize] = BASE_G;
    t[b'T' as usize] = BASE_T;
    t[b't' as usize] = BASE_T;
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

/// Classify the context of a monitored **top-strand** C, given the reference
/// base code immediately 3' of it (`ref[i+1]`). Returns `None` when that
/// neighbor is not a plain A/C/G/T (chrom end or ambiguity).
#[inline]
#[must_use]
pub(crate) fn top_context(next_code: u8) -> Option<Context> {
    match next_code {
        BASE_A => Some(Context::CpA),
        BASE_C => Some(Context::CpC),
        BASE_G => Some(Context::CpG),
        BASE_T => Some(Context::CpT),
        _ => None,
    }
}

/// Classify the context of a monitored **bottom-strand** C (i.e. a reference
/// `G` at position `i`), given the reference base code immediately 5' of the G
/// on the top strand (`ref[i-1]`). The cytosine sits on the bottom strand
/// complementary to the G, and its 3' neighbor on the bottom strand is the
/// complement of `ref[i-1]`. Returns `None` when `ref[i-1]` is not a plain
/// A/C/G/T.
///
/// Worked through: `ref[i-1] = C` → bottom context CpG; `T` → CpA; `A` → CpT;
/// `G` → CpC. (Matches NEB's "previous base == C ⇒ CpG" rule for monitored Gs.)
#[inline]
#[must_use]
pub(crate) fn bottom_context(prev_code: u8) -> Option<Context> {
    match prev_code {
        BASE_C => Some(Context::CpG),
        BASE_T => Some(Context::CpA),
        BASE_A => Some(Context::CpT),
        BASE_G => Some(Context::CpC),
        _ => None,
    }
}

// ── Reference access trait + encodings ──────────────────────────────────────

/// Per-contig reference accessor returning 4-bit BAM base codes by position.
///
/// Implementors back different in-memory layouts (1 byte/base, 4-bit packed,
/// 2-bit packed). The tally hot path is generic over this trait and
/// monomorphized per encoding, so there is no per-access branch on the layout.
pub(crate) trait RefCodes {
    /// Number of bases in the contig.
    fn len(&self) -> usize;
    /// The 4-bit BAM base code at `pos` (`pos < len`).
    fn code(&self, pos: usize) -> u8;
    /// Whether the base at `pos` equals the 4-bit BAM `code`. This is the
    /// reference-scan hot-path predicate; packed encodings override it to compare
    /// in their native representation and avoid decoding to the 4-bit code per
    /// base. `code` is a fixed monitored base (C=2 or G=4) across a scan.
    #[inline]
    fn monitors(&self, pos: usize, code: u8) -> bool {
        self.code(pos) == code
    }
    /// Top-strand context implied by the 3' neighbor base at `pos` (`ref[i+1]`).
    /// Packed encodings override to decode the neighbor in their native space.
    #[inline]
    fn ctx_top(&self, pos: usize) -> Option<Context> {
        top_context(self.code(pos))
    }
    /// Bottom-strand context implied by the 5' neighbor base at `pos`
    /// (`ref[i-1]`, complemented). Packed encodings override natively.
    #[inline]
    fn ctx_bottom(&self, pos: usize) -> Option<Context> {
        bottom_context(self.code(pos))
    }
}

/// One byte per base (the codes are stored directly). ~1 byte/base.
#[derive(Clone, Copy)]
pub(crate) struct ByteCodes<'a>(pub(crate) &'a [u8]);
impl RefCodes for ByteCodes<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.0.len()
    }
    #[inline]
    fn code(&self, pos: usize) -> u8 {
        self.0[pos]
    }
}

/// Two bases per byte (high nibble = even position, matching the BAM SEQ
/// convention). ~0.5 byte/base. Preserves the full 4-bit code space (incl. N).
#[derive(Clone, Copy)]
pub(crate) struct NibbleCodes<'a> {
    data: &'a [u8],
    len: usize,
}
impl RefCodes for NibbleCodes<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.len
    }
    #[inline]
    fn code(&self, pos: usize) -> u8 {
        let byte = self.data[pos >> 1];
        if pos & 1 == 0 { byte >> 4 } else { byte & 0x0F }
    }
}

/// Four bases per byte. ~0.25 byte/base. Only A/C/G/T are representable, so all
/// non-ACGT bases (N, IUPAC ambiguity) are folded to A at load. This preserves
/// monitored-site detection exactly (A is never the monitored C/G) but means a
/// monitored C/G adjacent to a former-N gets a concrete context instead of
/// being skipped (immaterial in practice — these sit in assembly gaps).
#[derive(Clone, Copy)]
pub(crate) struct TwoBitCodes<'a> {
    data: &'a [u8],
    len: usize,
}
impl RefCodes for TwoBitCodes<'_> {
    #[inline]
    fn len(&self) -> usize {
        self.len
    }
    #[inline]
    fn code(&self, pos: usize) -> u8 {
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        // 2-bit value → 4-bit BAM code without a table lookup: the codes for
        // A/C/G/T are 1/2/4/8 = 1 << (0/1/2/3).
        1u8 << val
    }
    #[inline]
    fn monitors(&self, pos: usize, code: u8) -> bool {
        // Compare in 2-bit space: `code == 1 << val`, so `val == log2(code)`.
        // This drops the per-base `1 << val` shift from the hot reference scan.
        // `code.trailing_zeros()` is loop-invariant (the monitored base), so it
        // hoists out of the caller's scan loop.
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        val == code.trailing_zeros() as u8
    }
    #[inline]
    fn ctx_top(&self, pos: usize) -> Option<Context> {
        // 2-bit value 0/1/2/3 (A/C/G/T) indexes the contexts CpA/CpC/CpG/CpT
        // directly, skipping the `1 << val` decode and the `match`. There is no N
        // in 2-bit (folded to A), so the context is always defined.
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        Some(Context::ALL[val as usize])
    }
    #[inline]
    fn ctx_bottom(&self, pos: usize) -> Option<Context> {
        // Bottom strand uses the complement of the neighbor: A↔T, C↔G, which on
        // the 2-bit value is `3 - val` (A0↔T3, C1↔G2).
        let val = (self.data[pos >> 2] >> ((pos & 3) * 2)) & 0x3;
        Some(Context::ALL[3 - val as usize])
    }
}

/// In-memory reference encoding, selected via `--ref-encoding`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RefEncoding {
    /// 1 byte per base (fastest, largest).
    Bytes,
    /// 4-bit packed, 2 bases per byte (~½ the memory).
    Nibble,
    /// 2-bit packed, 4 bases per byte (~¼ the memory; non-ACGT → A).
    TwoBit,
}

/// A packed contig: the packed bytes plus the base count (which the packing
/// can't recover from the byte length alone). Held inside [`Reference`]'s
/// packed variants; fields are private (access is via the typed accessors).
pub(crate) struct PackedContig {
    data: Vec<u8>,
    len: usize,
}

/// Fold a 4-bit code to its 2-bit value; non-ACGT (incl. N) → A (0).
#[inline]
fn nibble_to_2bit(code: u8) -> u8 {
    match code {
        BASE_C => 1,
        BASE_G => 2,
        BASE_T => 3,
        _ => 0, // A and anything non-ACGT
    }
}

fn pack_nibble(codes: &[u8]) -> PackedContig {
    let mut data = vec![0u8; codes.len().div_ceil(2)];
    for (i, &c) in codes.iter().enumerate() {
        if i & 1 == 0 {
            data[i >> 1] = c << 4;
        } else {
            data[i >> 1] |= c & 0x0F;
        }
    }
    PackedContig { data, len: codes.len() }
}

fn pack_twobit(codes: &[u8]) -> PackedContig {
    let mut data = vec![0u8; codes.len().div_ceil(4)];
    for (i, &c) in codes.iter().enumerate() {
        data[i >> 2] |= nibble_to_2bit(c) << ((i & 3) * 2);
    }
    PackedContig { data, len: codes.len() }
}

// ── Reference ───────────────────────────────────────────────────────────────

/// The full reference: every input `@SQ` contig loaded in the chosen encoding,
/// indexed by BAM `tid` (the 0-based position of the `@SQ` line).
pub(crate) enum Reference {
    /// 1 byte per base.
    Bytes(Vec<Vec<u8>>),
    /// 4-bit packed.
    Nibble(Vec<PackedContig>),
    /// 2-bit packed.
    TwoBit(Vec<PackedContig>),
}

impl Reference {
    /// Load every contig named by `header`'s `@SQ` lines from the indexed FASTA
    /// at `path`, in BAM tid order, packed per `encoding`.
    ///
    /// Semantics:
    /// * The FASTA may be a **superset** of the BAM's contigs — extra FASTA
    ///   contigs are ignored.
    /// * A `@SQ` with no corresponding FASTA entry is a fatal error.
    /// * A `@SQ` whose length disagrees with the FASTA `.fai` length is a fatal
    ///   error (the BAM was aligned against a different reference).
    ///
    /// Each contig is packed immediately after it is read, so peak load memory
    /// is one contig's worth of bytes plus the (smaller) accumulating packed
    /// store — not the full unpacked genome.
    ///
    /// # Errors
    /// Returns an error if the FASTA / `.fai` cannot be opened, a contig is
    /// missing, a length mismatches, or a contig cannot be read.
    pub(crate) fn load(
        path: &std::path::Path,
        header: &Header,
        encoding: RefEncoding,
    ) -> Result<Self> {
        let mut reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(path).with_context(
                || format!("opening indexed FASTA {} (is there a .fai?)", path.display()),
            )?;

        // Name → length from the .fai index, for the superset/length cross-check.
        let fai_lengths: HashMap<String, u64> = reader
            .index()
            .as_ref()
            .iter()
            .map(|rec| (String::from_utf8_lossy(rec.name().as_ref()).into_owned(), rec.length()))
            .collect();

        let n = header.reference_sequences().len();
        let mut bytes_v: Vec<Vec<u8>> =
            if encoding == RefEncoding::Bytes { Vec::with_capacity(n) } else { Vec::new() };
        let mut packed_v: Vec<PackedContig> =
            if encoding == RefEncoding::Bytes { Vec::new() } else { Vec::with_capacity(n) };

        for (name, map) in header.reference_sequences() {
            let name_str = std::str::from_utf8(name.as_ref())
                .map_err(|_| anyhow!("BAM @SQ name is not valid UTF-8"))?;
            let bam_len = usize::from(map.length()) as u64;

            match fai_lengths.get(name_str) {
                None => bail!(
                    "BAM contig '{name_str}' is not present in the reference FASTA index \
                     ({}). Every @SQ contig must exist in the reference.",
                    path.display()
                ),
                Some(&fai_len) if fai_len != bam_len => bail!(
                    "Length mismatch for contig '{name_str}': BAM @SQ says {bam_len} bp but the \
                     reference FASTA says {fai_len} bp. The BAM was aligned against a different \
                     reference."
                ),
                Some(_) => {}
            }

            // Seek to the contig's sequence start via the .fai byte offset and
            // stream it directly, avoiding the intermediate `Record` that
            // `query()` allocates (riker's trick — halves peak load memory).
            let region = Region::new(name_str, ..);
            let offset = reader.index().query(&region).with_context(|| {
                format!("looking up contig '{name_str}' in {}.fai", path.display())
            })?;
            reader
                .get_mut()
                .seek(SeekFrom::Start(offset))
                .with_context(|| format!("seeking to contig '{name_str}' in {}", path.display()))?;
            let mut raw = Vec::with_capacity(bam_len as usize);
            reader
                .read_sequence(&mut raw)
                .with_context(|| format!("reading contig '{name_str}' from {}", path.display()))?;
            let codes: Vec<u8> = raw.iter().map(|&b| REF_ASCII_TO_CODE[b as usize]).collect();
            match encoding {
                RefEncoding::Bytes => bytes_v.push(codes),
                RefEncoding::Nibble => packed_v.push(pack_nibble(&codes)),
                RefEncoding::TwoBit => packed_v.push(pack_twobit(&codes)),
            }
        }

        Ok(match encoding {
            RefEncoding::Bytes => Reference::Bytes(bytes_v),
            RefEncoding::Nibble => Reference::Nibble(packed_v),
            RefEncoding::TwoBit => Reference::TwoBit(packed_v),
        })
    }

    /// The active encoding.
    #[must_use]
    pub(crate) fn encoding(&self) -> RefEncoding {
        match self {
            Reference::Bytes(_) => RefEncoding::Bytes,
            Reference::Nibble(_) => RefEncoding::Nibble,
            Reference::TwoBit(_) => RefEncoding::TwoBit,
        }
    }

    /// Byte-encoded contig accessor for `tid` (only when encoding is `Bytes`).
    #[inline]
    #[must_use]
    pub(crate) fn byte_codes(&self, tid: i32) -> Option<ByteCodes<'_>> {
        match self {
            Reference::Bytes(v) if tid >= 0 => v.get(tid as usize).map(|c| ByteCodes(c)),
            _ => None,
        }
    }

    /// Nibble-packed contig accessor for `tid` (only when encoding is `Nibble`).
    #[inline]
    #[must_use]
    pub(crate) fn nibble_codes(&self, tid: i32) -> Option<NibbleCodes<'_>> {
        match self {
            Reference::Nibble(v) if tid >= 0 => {
                v.get(tid as usize).map(|c| NibbleCodes { data: &c.data, len: c.len })
            }
            _ => None,
        }
    }

    /// Two-bit-packed contig accessor for `tid` (only when encoding is `TwoBit`).
    #[inline]
    #[must_use]
    pub(crate) fn twobit_codes(&self, tid: i32) -> Option<TwoBitCodes<'_>> {
        match self {
            Reference::TwoBit(v) if tid >= 0 => {
                v.get(tid as usize).map(|c| TwoBitCodes { data: &c.data, len: c.len })
            }
            _ => None,
        }
    }

    /// Construct directly from already-encoded (byte) contigs (test helper).
    #[cfg(test)]
    #[must_use]
    pub(crate) fn from_encoded_contigs(contigs: Vec<Vec<u8>>) -> Self {
        Reference::Bytes(contigs)
    }
}

/// Encode an ASCII reference base to its 4-bit BAM code. Test-only: the load
/// path indexes [`REF_ASCII_TO_CODE`] directly; this is the named entry the
/// unit tests build encoded contigs with.
#[cfg(test)]
#[must_use]
pub(crate) fn encode_ref_base(ascii: u8) -> u8 {
    REF_ASCII_TO_CODE[ascii as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn top_context_buckets_by_next_base() {
        assert_eq!(top_context(BASE_A), Some(Context::CpA));
        assert_eq!(top_context(BASE_C), Some(Context::CpC));
        assert_eq!(top_context(BASE_G), Some(Context::CpG));
        assert_eq!(top_context(BASE_T), Some(Context::CpT));
        assert_eq!(top_context(BASE_N), None);
    }

    #[test]
    fn bottom_context_uses_complement_of_prev_base() {
        // ref[i-1]=C ⇒ CpG (the canonical CpG dinucleotide spanning i-1..i).
        assert_eq!(bottom_context(BASE_C), Some(Context::CpG));
        // ref[i-1]=T ⇒ complement A ⇒ CpA.
        assert_eq!(bottom_context(BASE_T), Some(Context::CpA));
        // ref[i-1]=A ⇒ complement T ⇒ CpT.
        assert_eq!(bottom_context(BASE_A), Some(Context::CpT));
        // ref[i-1]=G ⇒ complement C ⇒ CpC.
        assert_eq!(bottom_context(BASE_G), Some(Context::CpC));
        assert_eq!(bottom_context(BASE_N), None);
    }

    #[test]
    fn nibble_packing_round_trips_all_codes() {
        // Odd length to exercise the trailing half-byte.
        let codes: Vec<u8> = "CAGTNCCGTA".bytes().map(encode_ref_base).collect();
        let packed = pack_nibble(&codes);
        let view = NibbleCodes { data: &packed.data, len: packed.len };
        assert_eq!(view.len(), codes.len());
        for (i, &c) in codes.iter().enumerate() {
            assert_eq!(view.code(i), c, "nibble mismatch at {i}");
        }
        assert_eq!(packed.data.len(), codes.len().div_ceil(2));
    }

    #[test]
    fn twobit_packing_round_trips_acgt_and_folds_n_to_a() {
        let codes: Vec<u8> = "CAGTNCCGTA".bytes().map(encode_ref_base).collect();
        let packed = pack_twobit(&codes);
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
        let packed = pack_twobit(&codes);
        let view = TwoBitCodes { data: &packed.data, len: packed.len };
        for i in 0..view.len() {
            for &code in &[BASE_C, BASE_G] {
                assert_eq!(view.monitors(i, code), view.code(i) == code, "pos {i} code {code}");
            }
        }
    }

    #[test]
    fn context_index_is_stable() {
        assert_eq!(Context::CpA.index(), 0);
        assert_eq!(Context::CpC.index(), 1);
        assert_eq!(Context::CpG.index(), 2);
        assert_eq!(Context::CpT.index(), 3);
    }
}
