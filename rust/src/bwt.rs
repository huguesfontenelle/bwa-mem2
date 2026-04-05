/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of BWT / FM-index data structures from bwt.h and FMI_search.h.

Authors: Rust port based on original C/C++ by
         Sanchit Misra <sanchit.misra@intel.com>;
         Vasimuddin Md <vasimuddin.md@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

use anyhow::{bail, Context, Result};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

// ---------------------------------------------------------------------------
// Constants mirrored from FMI_search.h / macro.h
// ---------------------------------------------------------------------------

/// Number of BWT positions covered by one checkpoint entry.
pub const CP_BLOCK_SIZE: i64 = 64;

/// log2(CP_BLOCK_SIZE) — number of bits to shift right for the checkpoint index.
const CP_SHIFT: i64 = 6;

/// Mask to extract the within-block offset from a BWT position.
const CP_MASK: i64 = CP_BLOCK_SIZE - 1;

/// SA sampling interval (2^SA_COMPX = 2^3 = 8).
const SA_COMPX: i64 = 3;

/// Mask for checking whether a BWT position is a sampled SA entry.
const SA_COMPX_MASK: i64 = (1 << SA_COMPX) - 1; // 0x7

// ---------------------------------------------------------------------------
// Checkpoint occurrence entry
// ---------------------------------------------------------------------------

/// One checkpoint entry, covering 64 BWT positions.
///
/// Matches `CP_OCC` / `checkpoint_occ_scalar` in `FMI_search.h`.
#[derive(Debug, Clone, Default)]
pub struct CpOcc {
    /// Cumulative occurrence counts for A (0), C (1), G (2), T (3) at the
    /// start of this 64-base block.
    pub cp_count: [i64; 4],
    /// One-hot bit-vectors for each nucleotide in this 64-base block.
    /// Bit 63 (MSB) represents the first position in the block; bit 0 (LSB)
    /// represents the 64th position.  Only one of the four vectors can have
    /// a `1` at any given bit position (or all-zero for the sentinel `$`).
    pub one_hot_bwt_str: [u64; 4],
}

// ---------------------------------------------------------------------------
// FM-index
// ---------------------------------------------------------------------------

/// FM-index built from the 2-bit-packed BWT with 64-position checkpoints.
///
/// Matches the in-memory layout of `FMI_search` in `FMI_search.h`.
#[derive(Debug, Clone)]
pub struct FmIndex {
    /// Total length of the (doubled) reference sequence stored in the index.
    pub reference_seq_len: i64,
    /// BWT position of the sentinel character `$`.
    pub sentinel_index: i64,
    /// Checkpoint array — one entry per 64 BWT positions.
    pub cp_occ: Vec<CpOcc>,
    /// Lower 32 bits of the (compressed) suffix array.
    /// Entry `i` corresponds to BWT position `i * sa_intv`.
    pub sa_ls_word: Vec<u32>,
    /// Upper 8 bits of the (compressed) suffix array (signed, sign-extended
    /// to i64 when combined with `sa_ls_word`).
    pub sa_ms_byte: Vec<i8>,
    /// Prefix-sum table: `count[c]` is the number of BWT characters that are
    /// strictly less than `c` (after adding 1 for the sentinel).
    /// Indices: 0 = A, 1 = C, 2 = G, 3 = T, 4 = total (= seq_len + 1).
    pub count: [i64; 5],
    /// Suffix-array sampling interval (default 8 with `SA_COMPX = 3`).
    pub sa_intv: i64,
}

impl FmIndex {
    // -----------------------------------------------------------------------
    // Binary I/O
    // -----------------------------------------------------------------------

    /// Load the FM-index from `{prefix}.bwt.2bit.64`.
    ///
    /// File layout (all values little-endian):
    /// ```text
    /// i64          : reference_seq_len  (already incremented by 1 for the sentinel)
    /// i64 × 5      : raw count[0..5]    (before +1 adjustment)
    /// CpOcc × N    : checkpoint array   (N = (reference_seq_len >> 6) + 1)
    ///                  each CpOcc: i64×4 cp_count + u64×4 one_hot_bwt_str
    /// i8  × M      : sa_ms_byte         (M = (reference_seq_len >> 3) + 1)
    /// u32 × M      : sa_ls_word
    /// i64          : sentinel_index
    /// ```
    pub fn load(prefix: &str) -> Result<Self> {
        let bwt_path = format!("{}.bwt.2bit.64", prefix);
        let file = File::open(&bwt_path)
            .with_context(|| format!("cannot open BWT file '{}'", bwt_path))?;
        let mut r = BufReader::new(file);

        // --- reference_seq_len ---
        let reference_seq_len = r
            .read_i64::<LittleEndian>()
            .context("reading reference_seq_len")?;
        if reference_seq_len <= 0 {
            bail!("invalid reference_seq_len: {}", reference_seq_len);
        }

        // --- count[0..5] (raw; we add 1 to each after reading) ---
        let mut count = [0i64; 5];
        for c in &mut count {
            *c = r.read_i64::<LittleEndian>().context("reading count")?;
        }
        // The C loader increments every count[] entry by 1 after reading.
        for c in &mut count {
            *c += 1;
        }

        // --- checkpoint array ---
        let n_cp = (reference_seq_len >> CP_SHIFT) + 1;
        let mut cp_occ: Vec<CpOcc> = Vec::with_capacity(n_cp as usize);
        for _ in 0..n_cp {
            let mut entry = CpOcc::default();
            for v in &mut entry.cp_count {
                *v = r.read_i64::<LittleEndian>().context("reading cp_count")?;
            }
            for v in &mut entry.one_hot_bwt_str {
                *v = r.read_u64::<LittleEndian>().context("reading one_hot_bwt_str")?;
            }
            cp_occ.push(entry);
        }

        // --- compressed suffix array ---
        // n_sa = (reference_seq_len >> SA_COMPX) + 1
        let n_sa = (reference_seq_len >> SA_COMPX) + 1;
        let mut sa_ms_byte: Vec<i8> = Vec::with_capacity(n_sa as usize);
        for _ in 0..n_sa {
            sa_ms_byte.push(r.read_i8().context("reading sa_ms_byte")?);
        }
        let mut sa_ls_word: Vec<u32> = Vec::with_capacity(n_sa as usize);
        for _ in 0..n_sa {
            sa_ls_word.push(r.read_u32::<LittleEndian>().context("reading sa_ls_word")?);
        }

        // --- sentinel_index ---
        let sentinel_index = r
            .read_i64::<LittleEndian>()
            .context("reading sentinel_index")?;

        Ok(FmIndex {
            reference_seq_len,
            sentinel_index,
            cp_occ,
            sa_ls_word,
            sa_ms_byte,
            count,
            sa_intv: 1 << SA_COMPX,
        })
    }

    /// Save the FM-index to `{prefix}.bwt.2bit.64` in the same binary format
    /// that [`load`] expects.
    pub fn save(&self, prefix: &str) -> Result<()> {
        let bwt_path = format!("{}.bwt.2bit.64", prefix);
        let file = File::create(&bwt_path)
            .with_context(|| format!("cannot create BWT file '{}'", bwt_path))?;
        let mut w = BufWriter::new(file);

        // reference_seq_len
        w.write_i64::<LittleEndian>(self.reference_seq_len)
            .context("writing reference_seq_len")?;

        // count[0..5]: subtract 1 to reverse the +1 applied on load
        for &c in &self.count {
            w.write_i64::<LittleEndian>(c - 1)
                .context("writing count")?;
        }

        // checkpoint array
        for entry in &self.cp_occ {
            for &v in &entry.cp_count {
                w.write_i64::<LittleEndian>(v).context("writing cp_count")?;
            }
            for &v in &entry.one_hot_bwt_str {
                w.write_u64::<LittleEndian>(v)
                    .context("writing one_hot_bwt_str")?;
            }
        }

        // compressed suffix array
        for &b in &self.sa_ms_byte {
            w.write_i8(b).context("writing sa_ms_byte")?;
        }
        for &lw in &self.sa_ls_word {
            w.write_u32::<LittleEndian>(lw).context("writing sa_ls_word")?;
        }

        // sentinel_index
        w.write_i64::<LittleEndian>(self.sentinel_index)
            .context("writing sentinel_index")?;

        w.flush().context("flushing BWT file")?;
        Ok(())
    }

    // -----------------------------------------------------------------------
    // Core FM-index operations
    // -----------------------------------------------------------------------

    /// Count occurrences of character `c` in BWT[0..k] (exclusive upper bound).
    ///
    /// `c` is 0 for A, 1 for C, 2 for G, 3 for T.
    ///
    /// Implements the `GET_OCC` macro from `FMI_search.h`:
    /// ```text
    /// occ_id = k >> 6
    /// y      = k & 63          (within-block offset)
    /// occ    = cp_count[c]     (cumulative count at block start)
    ///        + popcount(one_hot_bwt_str[c] & mask_for_top_y_bits)
    /// ```
    pub fn occ(&self, k: i64, c: u8) -> i64 {
        debug_assert!(c < 4, "character index must be 0–3");
        if k < 0 {
            return 0;
        }

        let occ_id = (k >> CP_SHIFT) as usize;
        let y = (k & CP_MASK) as u32; // number of bits to count within the block

        let entry = &self.cp_occ[occ_id];
        let mut occ = entry.cp_count[c as usize];

        // The one-hot bit vector is stored MSB-first: bit 63 = position 0 of
        // the block, bit 62 = position 1, …, bit (63-y+1) = position y-1.
        //
        // `one_hot_mask_array[y]` from the C code has the top `y` bits set.
        // We reconstruct that mask here to avoid storing it separately.
        if y > 0 {
            // Shift right by (64 - y) to keep only the top y bits, then count.
            let bits = entry.one_hot_bwt_str[c as usize];
            let mask = if y >= 64 {
                !0u64
            } else {
                // Keep the top y bits: e.g. y=1 → 0x80..0, y=2 → 0xC0..0
                !0u64 << (64 - y)
            };
            occ += (bits & mask).count_ones() as i64;
        }

        occ
    }

    /// Return the BWT character at position `k` (0–3 for ACGT, 4 for sentinel).
    ///
    /// Uses the one-hot bit-vectors in the checkpoint that covers position `k`.
    pub fn bwt_char(&self, k: i64) -> u8 {
        debug_assert!(k >= 0 && k < self.reference_seq_len);

        let occ_id = (k >> CP_SHIFT) as usize;
        // Bit position within the u64: MSB = first position in block.
        // bit_pos = 63 − (k % 64): ranges from 63 (k % 64 == 0) to 0 (k % 64 == 63).
        let bit_pos = (CP_BLOCK_SIZE - 1 - (k & CP_MASK)) as u32;
        let entry = &self.cp_occ[occ_id];

        for c in 0u8..4 {
            if (entry.one_hot_bwt_str[c as usize] >> bit_pos) & 1 == 1 {
                return c;
            }
        }
        4 // sentinel '$'
    }

    /// LF-mapping: given BWT position `k` with character `c`, return the SA
    /// position of the preceding suffix.
    ///
    /// `lf(k, c) = count[c] + occ(k, c)`  (1-based occ counting up to k exclusive).
    pub fn lf_mapping(&self, k: i64, c: u8) -> i64 {
        self.count[c as usize] + self.occ(k, c)
    }

    /// One step of backward search: extend the BWT interval `[k, l]` with
    /// character `c` on the left.
    ///
    /// Returns `(k', l')` such that all suffixes in BWT[k'..=l'] are prefixed
    /// by `c` followed by the pattern that gave `[k, l]`.  Returns `(1, 0)`
    /// (an empty interval) if no extension exists.
    ///
    /// Formula:
    /// ```text
    /// k' = count[c] + occ(k − 1, c) + 1    (or count[c] + 1 if k == 0)
    /// l' = count[c] + occ(l, c)
    /// ```
    pub fn backward_extend(&self, k: i64, l: i64, c: u8) -> (i64, i64) {
        let occ_k_minus1 = if k > 0 { self.occ(k - 1, c) } else { 0 };
        let occ_l = self.occ(l, c);

        let base = self.count[c as usize];
        let new_k = base + occ_k_minus1 + 1;
        let new_l = base + occ_l;

        if new_k > new_l {
            (1, 0) // empty interval
        } else {
            (new_k, new_l)
        }
    }

    /// Look up the suffix array value at BWT position `i`.
    ///
    /// If `i` is a sampled position (`i % sa_intv == 0`), the value is read
    /// directly from the compressed arrays.  Otherwise LF-mapping is applied
    /// repeatedly until a sampled position is reached, and the offset is added.
    pub fn sa_lookup(&self, i: i64) -> i64 {
        let mask = SA_COMPX_MASK;

        if i & mask == 0 {
            // Direct lookup: sampled position.
            return self.sa_entry_at(i);
        }

        // Walk forward via LF-mapping until we hit a sampled position.
        let mut pos = i;
        let mut offset: i64 = 0;

        loop {
            let c = self.bwt_char(pos);
            if c == 4 {
                // Sentinel character: SA value is 0 by definition.
                return offset;
            }

            // LF-mapping: GET_OCC gives occ(pos, c)
            let occ_id = (pos >> CP_SHIFT) as usize;
            let y = (pos & CP_MASK) as u32;
            let entry = &self.cp_occ[occ_id];
            let occ = {
                let base_occ = entry.cp_count[c as usize];
                if y > 0 {
                    let bits = entry.one_hot_bwt_str[c as usize];
                    // Count bits in positions 0..y of the block (top y bits of u64).
                    let mask64 = if y >= 64 { !0u64 } else { !0u64 << (64 - y) };
                    base_occ + (bits & mask64).count_ones() as i64
                } else {
                    base_occ
                }
            };

            pos = self.count[c as usize] + occ;
            offset += 1;

            if pos & mask == 0 {
                break;
            }
        }

        self.sa_entry_at(pos) + offset
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /// Read the compressed suffix-array entry stored at *sampled* BWT position
    /// `pos` (i.e. `pos` must satisfy `pos % sa_intv == 0`).
    #[inline]
    fn sa_entry_at(&self, pos: i64) -> i64 {
        let idx = (pos >> SA_COMPX) as usize;
        let ms = self.sa_ms_byte[idx] as i64;
        let ls = self.sa_ls_word[idx] as i64;
        (ms << 32) | (ls & 0xffff_ffff)
    }
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a tiny synthetic FmIndex over the 4-base reference "ACGT".
    ///
    /// The doubled reference (forward + reverse complement) is "ACGTTGCA",
    /// length 8.  With the sentinel, `reference_seq_len` as stored in the file
    /// would be 9.  We hard-code small values just to exercise the arithmetic
    /// in `occ` and `bwt_char`.
    fn tiny_index() -> FmIndex {
        // We construct the simplest possible self-consistent index:
        // reference = "AA" (length 2), reference_seq_len = 3 (includes sentinel).
        // BWT = "A A $" at positions 0,1,2 but we won't verify BWT correctness;
        // we just want to exercise the Rust logic paths.

        // One checkpoint covers positions 0–63; we only have 3 BWT positions.
        // cp_count = cumulative counts before this block = all zeros.
        // one_hot_bwt_str: bit 63 = pos 0, bit 62 = pos 1, bit 61 = pos 2.
        // Suppose BWT = ['A', 'A', '$']:
        //   one_hot[A] = bits 63 and 62 set = 0xC000_0000_0000_0000
        //   one_hot[C] = 0, one_hot[G] = 0, one_hot[T] = 0

        let mut entry = CpOcc::default();
        entry.one_hot_bwt_str[0] = 0xC000_0000_0000_0000u64; // A at positions 0 and 1

        // count (after +1 adjustment): [1, 3, 3, 3, 3]
        // (0 As before pos 0, so count[A] = 0+1 = 1; no C/G/T at all = 3 after +1)
        let count = [1i64, 3, 3, 3, 3];

        // Compressed SA: sample every 8 positions.  With seq_len=3 we have 1 entry.
        let sa_ls_word = vec![0u32]; // SA[0] = 0 (position of smallest suffix)
        let sa_ms_byte = vec![0i8];

        FmIndex {
            reference_seq_len: 3,
            sentinel_index: 2, // '$' is at BWT position 2
            cp_occ: vec![entry],
            sa_ls_word,
            sa_ms_byte,
            count,
            sa_intv: 8,
        }
    }

    #[test]
    fn occ_zero_before_block() {
        let idx = tiny_index();
        // No characters before position 0 → occ(k<0, *) = 0
        assert_eq!(idx.occ(-1, 0), 0);
        assert_eq!(idx.occ(-1, 1), 0);
    }

    #[test]
    fn occ_counts_a_correctly() {
        let idx = tiny_index();
        // BWT pos 0 = A, pos 1 = A, pos 2 = '$'
        // occ(k, A): exclusive upper bound = k
        // occ(0, A): top 0 bits → 0
        assert_eq!(idx.occ(0, 0), 0);
        // occ(1, A): top 1 bit of 0xC0..0 = 1 bit → 1
        assert_eq!(idx.occ(1, 0), 1);
        // occ(2, A): top 2 bits of 0xC0..0 = 2 bits → 2
        assert_eq!(idx.occ(2, 0), 2);
    }

    #[test]
    fn occ_returns_zero_for_missing_chars() {
        let idx = tiny_index();
        // No C, G, T in our synthetic BWT.
        assert_eq!(idx.occ(2, 1), 0); // C
        assert_eq!(idx.occ(2, 2), 0); // G
        assert_eq!(idx.occ(2, 3), 0); // T
    }

    #[test]
    fn bwt_char_reads_correct_positions() {
        let idx = tiny_index();
        // Positions 0 and 1 are 'A' (0), position 2 is sentinel (4).
        assert_eq!(idx.bwt_char(0), 0); // A
        assert_eq!(idx.bwt_char(1), 0); // A
        assert_eq!(idx.bwt_char(2), 4); // sentinel
    }

    #[test]
    fn lf_mapping_basic() {
        let idx = tiny_index();
        // lf(k, A) = count[A] + occ(k, A) = 1 + occ(k, A)
        assert_eq!(idx.lf_mapping(0, 0), 1 + 0);
        assert_eq!(idx.lf_mapping(1, 0), 1 + 1);
        assert_eq!(idx.lf_mapping(2, 0), 1 + 2);
    }

    #[test]
    fn backward_extend_empty_on_no_match() {
        let idx = tiny_index();
        // If we extend with 'C' (1), there are no C's, so interval should be empty.
        let (k2, l2) = idx.backward_extend(1, 2, 1);
        assert!(k2 > l2, "expected empty interval, got ({}, {})", k2, l2);
    }

    #[test]
    fn backward_extend_a_in_small_index() {
        let idx = tiny_index();
        // Start with the full BWT interval [0, 2] (covering all 3 positions).
        // Extend with 'A': count[A]=1, occ(−1,A)=0, occ(2,A)=2
        // → k' = 1 + 0 + 1 = 2, l' = 1 + 2 = 3  but l'=3 > reference_seq_len-1=2
        // The logic itself is correct; validity of the interval depends on the
        // actual BWT.  We just assert the formula is applied as specified.
        let (k2, l2) = idx.backward_extend(0, 2, 0);
        assert_eq!(k2, 2);
        assert_eq!(l2, 3);
    }

    #[test]
    fn sa_lookup_sampled_position() {
        let idx = tiny_index();
        // Position 0 is a sampled position (0 % 8 == 0).
        // sa_ms_byte[0] = 0, sa_ls_word[0] = 0 → SA value = 0.
        assert_eq!(idx.sa_lookup(0), 0);
    }

    #[test]
    fn cp_block_size_constant() {
        assert_eq!(CP_BLOCK_SIZE, 64);
    }

    #[test]
    fn sa_intv_matches_compx() {
        let idx = tiny_index();
        assert_eq!(idx.sa_intv, 1 << SA_COMPX);
    }

    // Verify that save + load round-trips preserve all fields.
    #[test]
    fn save_load_roundtrip() {
        use std::fs;

        let original = tiny_index();
        let prefix = "/tmp/bwt_test_roundtrip";

        original.save(prefix).expect("save failed");

        let loaded = FmIndex::load(prefix).expect("load failed");

        assert_eq!(loaded.reference_seq_len, original.reference_seq_len);
        assert_eq!(loaded.sentinel_index, original.sentinel_index);
        assert_eq!(loaded.sa_intv, original.sa_intv);
        assert_eq!(loaded.count, original.count);
        assert_eq!(loaded.sa_ls_word, original.sa_ls_word);
        assert_eq!(loaded.sa_ms_byte, original.sa_ms_byte);
        assert_eq!(loaded.cp_occ.len(), original.cp_occ.len());
        for (a, b) in loaded.cp_occ.iter().zip(original.cp_occ.iter()) {
            assert_eq!(a.cp_count, b.cp_count);
            assert_eq!(a.one_hot_bwt_str, b.one_hot_bwt_str);
        }

        // Clean up.
        let _ = fs::remove_file(format!("{}.bwt.2bit.64", prefix));
    }
}
