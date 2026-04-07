/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of FM-index construction from FMI_search.cpp and bwtindex.cpp.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

use crate::bwt::{FmIndex, CpOcc};
use crate::bntseq::BntSeq;
use anyhow::{Result, bail, Context};

// ---------------------------------------------------------------------------
// SA-IS: Suffix Array Induced Sorting
// Reference: Nong, Zhang, Chan (2009)
// ---------------------------------------------------------------------------

/// Classify each position as S-type (true) or L-type (false).
/// The last character (sentinel) is always S-type.
fn classify_sl(s: &[u8]) -> Vec<bool> {
    let n = s.len();
    let mut t = vec![false; n];
    if n == 0 {
        return t;
    }
    // The last character is the sentinel (smallest), always S-type.
    t[n - 1] = true;
    if n == 1 {
        return t;
    }
    // Scan right to left.
    for i in (0..n - 1).rev() {
        t[i] = if s[i] < s[i + 1] {
            true // S-type
        } else if s[i] > s[i + 1] {
            false // L-type
        } else {
            t[i + 1] // same type as next
        };
    }
    t
}

/// Returns true if position i is a Left-Most S (LMS) suffix.
#[inline]
fn is_lms(t: &[bool], i: usize) -> bool {
    i > 0 && t[i] && !t[i - 1]
}

/// Build bucket head (first index of each bucket) or tail (last index).
fn get_buckets(s: &[u8], alphabet_size: usize, end: bool) -> Vec<usize> {
    let mut bkt_sizes = vec![0usize; alphabet_size];
    for &c in s.iter() {
        bkt_sizes[c as usize] += 1;
    }
    let mut buckets = vec![0usize; alphabet_size];
    let mut sum = 0usize;
    for i in 0..alphabet_size {
        sum += bkt_sizes[i];
        buckets[i] = if end { sum - 1 } else { sum - bkt_sizes[i] };
    }
    buckets
}

/// Induced sorting: scan left to right, place L-type suffixes.
fn induced_sort_l(
    s: &[u8],
    sa: &mut Vec<i64>,
    t: &[bool],
    alphabet_size: usize,
) {
    let mut bkt = get_buckets(s, alphabet_size, false); // heads
    let n = s.len();
    for i in 0..n {
        let j = sa[i];
        if j <= 0 {
            continue;
        }
        let k = (j - 1) as usize;
        if !t[k] {
            sa[bkt[s[k] as usize]] = k as i64;
            bkt[s[k] as usize] += 1;
        }
    }
}

/// Induced sorting: scan right to left, place S-type suffixes.
fn induced_sort_s(
    s: &[u8],
    sa: &mut Vec<i64>,
    t: &[bool],
    alphabet_size: usize,
) {
    let mut bkt = get_buckets(s, alphabet_size, true); // tails
    let n = s.len();
    for i in (0..n).rev() {
        let j = sa[i];
        if j <= 0 {
            continue;
        }
        let k = (j - 1) as usize;
        if t[k] {
            sa[bkt[s[k] as usize]] = k as i64;
            if bkt[s[k] as usize] > 0 {
                bkt[s[k] as usize] -= 1;
            }
        }
    }
}

/// Core SA-IS recursion.  `alphabet_size` is the number of distinct characters.
fn sais_recursive(s: &[u8], alphabet_size: usize) -> Vec<i64> {
    let n = s.len();
    let mut sa = vec![-1i64; n];

    let t = classify_sl(s);

    // Collect LMS positions.
    let lms_positions: Vec<usize> = (1..n).filter(|&i| is_lms(&t, i)).collect();

    // Step 1: Place LMS suffixes into tail of their buckets.
    {
        let mut bkt = get_buckets(s, alphabet_size, true); // tails
        for &pos in lms_positions.iter().rev() {
            sa[bkt[s[pos] as usize]] = pos as i64;
            if bkt[s[pos] as usize] > 0 {
                bkt[s[pos] as usize] -= 1;
            }
        }
    }

    // Step 2: Induce L-type then S-type.
    induced_sort_l(s, &mut sa, &t, alphabet_size);
    induced_sort_s(s, &mut sa, &t, alphabet_size);

    // Step 3: Compact sorted LMS substrings and assign names.
    // Collect all LMS positions in sorted order.
    let sorted_lms: Vec<usize> = sa.iter()
        .filter(|&&x| x > 0 && is_lms(&t, x as usize))
        .map(|&x| x as usize)
        .collect();

    // Assign lexicographic names to LMS substrings.
    let mut name = 0i64;
    let mut prev: Option<usize> = None;
    let mut lms_names = vec![-1i64; n]; // indexed by original position

    for &pos in &sorted_lms {
        // Check if this LMS substring differs from the previous one.
        let mut differs = true;
        if let Some(p) = prev {
            // Compare character by character until both substrings end.
            let mut k = 0usize;
            loop {
                let same = (s[pos + k] == s[p + k]) && (t[pos + k] == t[p + k]);
                let end_pos = k > 0 && is_lms(&t, pos + k);
                let end_p = k > 0 && is_lms(&t, p + k);
                if !same || end_pos || end_p {
                    differs = !(same && end_pos && end_p);
                    break;
                }
                k += 1;
            }
        }
        if differs {
            name += 1;
            prev = Some(pos);
        }
        lms_names[pos] = name - 1;
    }

    let num_names = name as usize;

    // Build reduced string: only LMS positions, in order of occurrence.
    let reduced: Vec<u8> = lms_positions.iter()
        .map(|&pos| lms_names[pos] as u8)
        .collect();

    // Solve reduced problem.
    let reduced_sa: Vec<i64> = if num_names < reduced.len() {
        // Names are not all unique; recurse.
        sais_recursive(&reduced, num_names)
    } else {
        // All names unique; build SA directly by inverting names.
        let mut rsa = vec![0i64; reduced.len()];
        for (i, &name) in reduced.iter().enumerate() {
            rsa[name as usize] = i as i64;
        }
        rsa
    };

    // Step 4: Induce final SA from reduced SA.
    // Re-scatter LMS suffixes in the order given by reduced_sa.
    sa.iter_mut().for_each(|x| *x = -1);
    {
        let mut bkt = get_buckets(s, alphabet_size, true); // tails
        // Place in reverse order so earlier ones end up at higher indices first.
        for i in (0..reduced_sa.len()).rev() {
            let lms_pos = lms_positions[reduced_sa[i] as usize];
            sa[bkt[s[lms_pos] as usize]] = lms_pos as i64;
            if bkt[s[lms_pos] as usize] > 0 {
                bkt[s[lms_pos] as usize] -= 1;
            }
        }
    }

    induced_sort_l(s, &mut sa, &t, alphabet_size);
    induced_sort_s(s, &mut sa, &t, alphabet_size);

    sa
}

/// Build the suffix array of `s` using SA-IS.
///
/// The input alphabet for DNA is `{0,1,2,3}` plus sentinel `4`.
/// Returns `Vec<i64>` of length `s.len()`.
///
/// # Sentinel convention
///
/// SA-IS requires the sentinel to be **strictly smaller** than every other
/// character.  In our DNA+sentinel encoding the sentinel has value `4`
/// (the *largest* value), so we remap the alphabet before invoking the
/// recursion:
///
/// ```text
/// sentinel (max_val) → 0
/// every other character c → c + 1
/// ```
///
/// This preserves relative ordering among non-sentinel characters and
/// guarantees that the sentinel sorts first.  The suffix-array positions
/// returned are indices into the *original* (un-remapped) string, so callers
/// do not need to adjust them.
pub fn build_suffix_array(s: &[u8]) -> Vec<i64> {
    if s.is_empty() {
        return Vec::new();
    }
    let max_val = *s.iter().max().unwrap() as usize;
    // Remap: sentinel (max_val) → 0; every other value c → c + 1.
    let remapped: Vec<u8> = s.iter()
        .map(|&c| if c as usize == max_val { 0u8 } else { c + 1 })
        .collect();
    sais_recursive(&remapped, max_val + 1)
}

// ---------------------------------------------------------------------------
// BWT from SA
// ---------------------------------------------------------------------------

/// Construct the BWT from the suffix array.
///
/// `bwt[i] = s[sa[i] - 1]` if `sa[i] > 0`, else the sentinel value `4`.
pub fn build_bwt_from_sa(s: &[u8], sa: &[i64]) -> Vec<u8> {
    let n = sa.len();
    let mut bwt = vec![0u8; n];
    for i in 0..n {
        if sa[i] == 0 {
            bwt[i] = 4; // sentinel '$'
        } else {
            bwt[i] = s[(sa[i] - 1) as usize];
        }
    }
    bwt
}

// ---------------------------------------------------------------------------
// FM-index construction
// ---------------------------------------------------------------------------

/// Build the FM-index from a packed DNA sequence (values in `{0,1,2,3}`).
///
/// The sequence should NOT include the sentinel; it is appended internally.
pub fn build_fm_index(seq: &[u8], seq_len: i64) -> Result<FmIndex> {
    let n = seq_len as usize;
    if n == 0 {
        bail!("build_fm_index: empty sequence");
    }
    if n != seq.len() {
        bail!("build_fm_index: seq_len ({}) != seq.len() ({})", n, seq.len());
    }

    // ------------------------------------------------------------------
    // 1. Append sentinel (4) and build the suffix array.
    // ------------------------------------------------------------------
    let mut s_with_sentinel: Vec<u8> = Vec::with_capacity(n + 1);
    s_with_sentinel.extend_from_slice(seq);
    s_with_sentinel.push(4u8); // sentinel

    eprintln!("[build_fm_index] Building suffix array for {} bp sequence ...", n);
    let sa = build_suffix_array(&s_with_sentinel);
    let total_len = sa.len(); // n + 1

    // ------------------------------------------------------------------
    // 2. Find the sentinel index (position in BWT where SA[i] == 0).
    // ------------------------------------------------------------------
    let sentinel_index = sa.iter()
        .position(|&x| x == 0)
        .context("build_fm_index: sentinel not found in suffix array")? as i64;

    // ------------------------------------------------------------------
    // 3. Build BWT from SA.
    // ------------------------------------------------------------------
    eprintln!("[build_fm_index] Building BWT ...");
    let bwt = build_bwt_from_sa(&s_with_sentinel, &sa);

    // ------------------------------------------------------------------
    // 4. Compute count[] (cumulative nucleotide counts).
    //    count[c] = number of BWT characters strictly less than c
    //    (treating sentinel as character index -1, i.e., less than 0).
    //    Layout: count[0] = 1 (one sentinel before all A's).
    //            count[c+1] = count[c] + number of c in BWT, for c in 0..4
    // ------------------------------------------------------------------
    let mut raw_counts = [0i64; 5]; // raw_counts[c] = occurrences of c in BWT
    for &b in bwt.iter() {
        if b < 4 {
            raw_counts[b as usize] += 1;
        }
        // sentinel (b==4) is counted implicitly: there is exactly one.
    }

    let mut count = [0i64; 5];
    // count[0]: number of characters strictly less than A = 1 (one sentinel '$').
    count[0] = 1;
    // count[c+1] = count[c] + occurrences of character c in BWT.
    // After the loop: count[c] = number of BWT characters < c (including sentinel).
    for c in 0..4usize {
        count[c + 1] = count[c] + raw_counts[c];
    }

    eprintln!(
        "[build_fm_index] count = [{}, {}, {}, {}, {}]",
        count[0], count[1], count[2], count[3], count[4]
    );

    // ------------------------------------------------------------------
    // 5. Build checkpoints (one per CP_BLOCK_SIZE = 64 BWT positions).
    //    Each CpOcc stores:
    //      cp_count[c]          = cumulative count of c BEFORE this block
    //      one_hot_bwt_str[c]   = 64-bit one-hot vector for character c
    //                             across the 64 positions in this block.
    //
    //    one_hot_bwt_str is built MSB-first: bit 63 is bwt[block_start+0],
    //    bit 0 is bwt[block_start+63].  This matches the GET_OCC macro in
    //    FMI_search.h which masks with one_hot_mask_array[y] where mask[1]
    //    has bit 63 set, mask[2] has bits 63:62 set, etc.
    // ------------------------------------------------------------------
    const CP_BLOCK_SIZE: usize = 64;
    const CP_SHIFT: usize = 6; // log2(64)

    // Align total_len up to a CP_BLOCK_SIZE multiple for safe block reads.
    let total_len_aligned = ((total_len + CP_BLOCK_SIZE - 1) / CP_BLOCK_SIZE) * CP_BLOCK_SIZE;

    // Build a padded BWT array (padding with DUMMY_CHAR=6 beyond total_len).
    let mut bwt_padded = bwt.clone();
    bwt_padded.resize(total_len_aligned, 6u8);

    let cp_occ_size = (total_len >> CP_SHIFT) + 1;
    let mut cp_occ: Vec<CpOcc> = Vec::with_capacity(cp_occ_size);

    let mut cp_count = [0i64; 4]; // running counts of {A,C,G,T}

    for i in 0..total_len {
        if (i & (CP_BLOCK_SIZE - 1)) == 0 {
            // Start of a new block: record cumulative counts and build one-hot vectors.
            let block_start = i;
            let mut one_hot = [0u64; 4];

            for j in 0..CP_BLOCK_SIZE {
                // Shift existing bits left to make room for the new bit at LSB.
                for c in 0..4usize {
                    one_hot[c] <<= 1;
                }
                let ch = bwt_padded[block_start + j];
                if ch < 4 {
                    one_hot[ch as usize] |= 1;
                }
            }

            cp_occ.push(CpOcc {
                cp_count: [cp_count[0], cp_count[1], cp_count[2], cp_count[3]],
                one_hot_bwt_str: one_hot,
            });
        }

        // Update running counts.
        let ch = bwt_padded[i];
        if ch < 4 {
            cp_count[ch as usize] += 1;
        }
    }

    // ------------------------------------------------------------------
    // 6. Build compressed suffix array.
    //    SA sampling interval: sa_intv = 32.
    //    For every i where sa[i] % sa_intv == 0, store the SA value split
    //    into lower 32 bits (sa_ls_word) and upper byte (sa_ms_byte).
    //    Storage is indexed by i / sa_intv (NOT by sa[i] / sa_intv).
    //    This matches the C code: samples are taken at every sa_intv-th
    //    BWT position, not every sa_intv-th SA value.
    // ------------------------------------------------------------------
    const SA_INTV: i64 = 32;
    let sa_count = ((total_len as i64) / SA_INTV + 1) as usize;

    let mut sa_ls_word: Vec<u32> = vec![0u32; sa_count];
    let mut sa_ms_byte: Vec<i8> = vec![0i8; sa_count];

    let mut pos = 0usize;
    for i in 0..total_len {
        if (i as i64 % SA_INTV) == 0 {
            let val = sa[i];
            sa_ls_word[pos] = (val & 0xffff_ffff) as u32;
            sa_ms_byte[pos] = ((val >> 32) & 0xff) as i8;
            pos += 1;
        }
    }

    eprintln!("[build_fm_index] FM-index built successfully.");

    Ok(FmIndex {
        reference_seq_len: total_len as i64,
        sentinel_index,
        cp_occ,
        sa_ls_word,
        sa_ms_byte,
        count,
        sa_intv: SA_INTV,
    })
}

// ---------------------------------------------------------------------------
// Index building entry point
// ---------------------------------------------------------------------------

/// Build the BWA-MEM2 FM-index for a FASTA file.
///
/// Steps:
/// 1. Pack the FASTA into a binary sequence via `BntSeq::fasta_to_bntseq`.
/// 2. Build the FM-index from the packed sequence.
/// 3. Save the index to disk via `FmIndex::save`.
pub fn build_index(fasta: &str, prefix: &str) -> Result<()> {
    eprintln!("[build_index] Packing FASTA ...");
    let (bntseq, pac) = BntSeq::fasta_to_bntseq(fasta, prefix)
        .with_context(|| format!("Failed to pack FASTA '{}'", fasta))?;

    let seq_len = bntseq.l_pac;

    eprintln!("[build_index] Reference length: {} bp", seq_len);

    if seq_len == 0 {
        bail!("build_index: reference sequence has zero length");
    }

    // Unpack from 4-bases-per-byte to 1-base-per-byte (0..3)
    let seq: Vec<u8> = (0..seq_len).map(|i| {
        let byte_idx = (i >> 2) as usize;
        let shift = ((3 - (i & 3)) * 2) as u32;
        if byte_idx < pac.len() { (pac[byte_idx] >> shift) & 3 } else { 0 }
    }).collect();

    eprintln!("[build_index] Building FM-index ...");
    let fm_index = build_fm_index(&seq, seq_len)
        .context("Failed to build FM-index")?;

    eprintln!("[build_index] Saving FM-index to disk ...");
    fm_index.save(prefix)
        .with_context(|| format!("Failed to save FM-index with prefix '{}'", prefix))?;

    eprintln!("[build_index] Done.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Verify that the suffix array of a small string is correct.
    #[test]
    fn sa_is_small() {
        // DNA sequence "ACGT" + sentinel (4).
        // Encoded: [0, 1, 2, 3, 4]
        // Suffixes (0-indexed):
        //   0: ACGT$  -> SA entry 0
        //   1: CGT$   -> SA entry 1
        //   2: GT$    -> SA entry 2
        //   3: T$     -> SA entry 3
        //   4: $      -> SA entry 4 (smallest)
        // Sorted order: $ < A < C < G < T
        // SA = [4, 0, 1, 2, 3]
        let s: Vec<u8> = vec![0, 1, 2, 3, 4];
        let sa = build_suffix_array(&s);
        assert_eq!(sa, vec![4, 0, 1, 2, 3]);
    }

    /// Verify SA for a repeated pattern.
    #[test]
    fn sa_is_repeated() {
        // "AABB" + sentinel: [0,0,1,1,4]  (A=0, B=1, $=4)
        //
        // Treating $ as the SMALLEST character the sorted suffixes are:
        //   pos 4: [$]                 → name "$"
        //   pos 0: [A,A,B,B,$]        → "AABB$"
        //   pos 1: [A,B,B,$]          → "ABB$"
        //   pos 3: [B,$]              → "B$"    (B$ < BB$ because $ < B)
        //   pos 2: [B,B,$]            → "BB$"
        // Expected SA = [4, 0, 1, 3, 2]
        let s: Vec<u8> = vec![0, 0, 1, 1, 4];
        let sa = build_suffix_array(&s);
        assert_eq!(sa, vec![4, 0, 1, 3, 2],
            "SA mismatch for [0,0,1,1,4]: got {:?}", sa);

        // Also verify ordering with a sentinel-aware comparator.
        // Remap so that the sentinel (4) sorts as 0 and other chars shift up.
        let remap = |c: u8| if c == 4 { 0u8 } else { c + 1 };
        let remapped: Vec<u8> = s.iter().map(|&c| remap(c)).collect();
        let n = s.len();
        for i in 1..n {
            let a = sa[i - 1] as usize;
            let b = sa[i] as usize;
            assert!(
                remapped[a..] <= remapped[b..],
                "SA not sorted at index {}: sa[{}]={} vs sa[{}]={}",
                i, i - 1, a, i, b
            );
        }
    }

    /// BWT construction sanity check.
    #[test]
    fn bwt_from_sa_sentinel() {
        // "ACGT" + sentinel (4)
        let s: Vec<u8> = vec![0, 1, 2, 3, 4];
        let sa = build_suffix_array(&s);
        let bwt = build_bwt_from_sa(&s, &sa);
        // sa = [4, 0, 1, 2, 3]
        // bwt[0] = s[4-1] = s[3] = 3  (T)
        // bwt[1] = sentinel since sa[1]=0 -> 4
        // bwt[2] = s[0] = 0 (A)
        // bwt[3] = s[1] = 1 (C)
        // bwt[4] = s[2] = 2 (G)
        assert_eq!(bwt, vec![3, 4, 0, 1, 2]);
    }

    /// Count array must be monotonically non-decreasing.
    #[test]
    fn count_array_monotone() {
        // "AACGT" encoded: [0,0,1,2,3]
        let seq: Vec<u8> = vec![0, 0, 1, 2, 3];
        let fm = build_fm_index(&seq, 5).expect("build_fm_index failed");
        // count should be non-decreasing
        for i in 1..5 {
            assert!(
                fm.count[i] >= fm.count[i - 1],
                "count not monotone: count[{}]={} < count[{}]={}",
                i, fm.count[i], i-1, fm.count[i-1]
            );
        }
        // count[0] must equal 1 (one sentinel)
        assert_eq!(fm.count[0], 1);
    }

    /// Checkpoint count values at block boundaries should match a manual scan.
    #[test]
    fn cp_occ_counts_correct() {
        // Use a simple 65-character sequence so we get two checkpoints.
        let seq: Vec<u8> = (0..65u8).map(|i| i % 4).collect();
        let fm = build_fm_index(&seq, 65).expect("build_fm_index failed");
        // cp_occ[0].cp_count should all be zero (no characters before block 0).
        let first = &fm.cp_occ[0];
        for c in 0..4 {
            assert_eq!(first.cp_count[c], 0, "cp_count[{}] at block 0 should be 0", c);
        }
    }

    /// SA interval must be 32.
    #[test]
    fn sa_intv_is_32() {
        let seq: Vec<u8> = vec![0, 1, 2, 3, 0, 1, 2, 3];
        let fm = build_fm_index(&seq, 8).expect("build_fm_index failed");
        assert_eq!(fm.sa_intv, 32);
    }

    /// The sentinel index must be within range.
    #[test]
    fn sentinel_index_in_range() {
        let seq: Vec<u8> = vec![0, 1, 2, 3];
        let fm = build_fm_index(&seq, 4).expect("build_fm_index failed");
        assert!(
            fm.sentinel_index >= 0 && fm.sentinel_index < fm.reference_seq_len,
            "sentinel_index {} out of range [0, {})",
            fm.sentinel_index, fm.reference_seq_len
        );
    }
}
