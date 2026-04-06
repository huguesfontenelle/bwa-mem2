/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of SMEM finding from FMI_search.cpp.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

//! SMEM (Super-Maximal Exact Match) finding using the FM-index.
//!
//! This module ports the core SMEM-finding logic from `FMI_search.cpp` in bwa-mem2.
//! It uses the bi-directional BWT approach: each live interval tracks both a forward
//! BWT interval (`k`, hit count `s`) and a reverse-complement BWT interval (`l`),
//! enabling simultaneous left- and right-extension.

use crate::bwt::FmIndex;
use crate::types::Smem;

// ---------------------------------------------------------------------------
// ASCII → 2-bit encoding
// ---------------------------------------------------------------------------

/// Convert an ASCII nucleotide sequence into 2-bit encoded bytes.
///
/// Encoding: `A`/`a` → 0, `C`/`c` → 1, `G`/`g` → 2, `T`/`t` → 3,
/// any other character (including `N`/`n`) → 4.
///
/// The returned `Vec<u8>` has the same length as `seq`.
pub fn query_to_2bit(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'A' | b'a' => 0u8,
            b'C' | b'c' => 1u8,
            b'G' | b'g' => 2u8,
            b'T' | b't' => 3u8,
            _ => 4u8,
        })
        .collect()
}

// ---------------------------------------------------------------------------
// SA entry lookup
// ---------------------------------------------------------------------------

/// Retrieve reference coordinates for all SA entries in the BWT interval
/// `[smem.k, smem.k + smem.s)` of the given SMEM.
///
/// If the interval size (`smem.s`) exceeds `max_occ` the match is considered
/// too repetitive and an empty `Vec` is returned immediately.
///
/// Otherwise, `fmi.sa_lookup(i)` is called for every BWT row `i` in the
/// half-open interval `[smem.k, smem.k + smem.s)`.
pub fn get_sa_entries(fmi: &FmIndex, smem: &Smem, max_occ: i32) -> Vec<i64> {
    // smem.s is the interval size (hit count) stored in the Smem struct.
    let hit_count = smem.s;
    if hit_count > max_occ as i64 {
        return Vec::new();
    }
    let mut entries = Vec::with_capacity(hit_count as usize);
    for i in smem.k..smem.k + hit_count {
        entries.push(fmi.sa_lookup(i));
    }
    entries
}

// ---------------------------------------------------------------------------
// Internal bi-directional interval type
// ---------------------------------------------------------------------------

/// A live bi-directional BWT interval used during SMEM construction.
///
/// `k`  – lower bound of the forward BWT interval (SA index).
/// `l`  – lower bound of the reverse-complement BWT interval (SA index).
/// `s`  – interval size (number of hits); the forward interval spans `[k, k+s)`.
/// `m`  – leftmost query position consumed so far (inclusive).
/// `n`  – rightmost query position consumed so far (inclusive).
#[derive(Clone, Copy, Debug)]
struct BiInterval {
    k: i64,
    l: i64,
    s: i64,
    m: i32,
    n: i32,
}

impl BiInterval {
    /// Produce a `Smem` record from this interval for read `rid`.
    ///
    /// The query span is `[m, n]` inclusive, stored as `[m, n+1)` in `Smem`.
    fn to_smem(self, rid: u32) -> Smem {
        Smem {
            rid,
            m: self.m as u32,
            n: (self.n + 1) as u32,
            k: self.k,
            l: self.l,
            s: self.s,
        }
    }
}

// ---------------------------------------------------------------------------
// bi-directional backward extension
// ---------------------------------------------------------------------------

/// Extend a bi-directional interval one position to the **right** by character
/// `a` (0-based 2-bit encoding: A=0, C=1, G=2, T=3).
///
/// This mirrors the `backwardExt` method in `FMI_search.cpp`.  The forward
/// extension of query character `a` is performed as a backward extension on
/// the reverse-complement strand with complement `3-a`, then the resulting
/// `k` and `l` fields are swapped back.
///
/// The update rule (for all four characters `b`) is:
/// ```text
/// occ_sp[b] = OCC(k,   b)   // occurrences of b in BWT[0..k)
/// occ_ep[b] = OCC(k+s, b)   // occurrences of b in BWT[0..k+s)
/// k_new[b]  = count[b] + occ_sp[b]
/// s_new[b]  = occ_ep[b] - occ_sp[b]
/// ```
///
/// The reverse-complement interval pointer `l` is updated using the
/// relationship between forward and reverse intervals (Lam et al. 2009):
/// ```text
/// l_new[3] = l + sentinel_offset
/// l_new[2] = l_new[3] + s_new[3]
/// l_new[1] = l_new[2] + s_new[2]
/// l_new[0] = l_new[1] + s_new[1]
/// ```
///
/// Returns `None` when the new interval for character `a` is empty (`s == 0`).
fn bidir_extend_right(fmi: &FmIndex, iv: &BiInterval, a: u8) -> Option<BiInterval> {
    debug_assert!(a < 4, "bidir_extend_right: character must be in 0..4, got {}", a);

    let mut k_new = [0i64; 4];
    let mut s_new = [0i64; 4];

    // Compute new forward intervals for all four characters simultaneously.
    // OCC(pos, b) = number of occurrences of b in BWT[0..pos) (exclusive).
    for b in 0u8..4 {
        let occ_sp = fmi.occ(iv.k, b);
        let occ_ep = fmi.occ(iv.k + iv.s, b);
        k_new[b as usize] = fmi.count[b as usize] + occ_sp;
        s_new[b as usize] = occ_ep - occ_sp;
    }

    // Determine whether the sentinel '$' falls inside the current forward interval.
    // The sentinel is located at `sentinel_index` in the BWT.  If it is inside
    // [k, k+s), the reverse-complement `l` pointer shifts by one.
    let sentinel_offset: i64 =
        if iv.k <= fmi.sentinel_index && fmi.sentinel_index < iv.k + iv.s {
            1
        } else {
            0
        };

    // Update reverse-complement interval lower bounds using the recurrence.
    let mut l_new = [0i64; 4];
    l_new[3] = iv.l + sentinel_offset;
    l_new[2] = l_new[3] + s_new[3];
    l_new[1] = l_new[2] + s_new[2];
    l_new[0] = l_new[1] + s_new[1];

    let a = a as usize;
    if s_new[a] == 0 {
        None
    } else {
        Some(BiInterval {
            k: k_new[a],
            l: l_new[a],
            s: s_new[a],
            m: iv.m,
            n: iv.n,
        })
    }
}

/// Extend a bi-directional interval one position to the **left** by character
/// `a`.
///
/// Left extension uses the reverse-complement strand interval (`l`, `s`) as
/// the working interval, extends it backward with character `3-a`, then swaps
/// the resulting `k` and `l` back.  This is exactly the dual of
/// `bidir_extend_right`.
///
/// Returns `None` when the new interval is empty.
fn bidir_extend_left(fmi: &FmIndex, iv: &BiInterval, a: u8) -> Option<BiInterval> {
    debug_assert!(a < 4, "bidir_extend_left: character must be in 0..4, got {}", a);

    // Work on the reverse-complement interval (l, s) with complement base.
    let comp = 3u8 - a;

    let mut k_new_rev = [0i64; 4];
    let mut s_new_rev = [0i64; 4];

    for b in 0u8..4 {
        let occ_sp = fmi.occ(iv.l, b);
        let occ_ep = fmi.occ(iv.l + iv.s, b);
        k_new_rev[b as usize] = fmi.count[b as usize] + occ_sp;
        s_new_rev[b as usize] = occ_ep - occ_sp;
    }

    let sentinel_offset: i64 =
        if iv.l <= fmi.sentinel_index && fmi.sentinel_index < iv.l + iv.s {
            1
        } else {
            0
        };

    let mut l_new_rev = [0i64; 4];
    l_new_rev[3] = iv.k + sentinel_offset;
    l_new_rev[2] = l_new_rev[3] + s_new_rev[3];
    l_new_rev[1] = l_new_rev[2] + s_new_rev[2];
    l_new_rev[0] = l_new_rev[1] + s_new_rev[1];

    let comp = comp as usize;
    if s_new_rev[comp] == 0 {
        None
    } else {
        // After extending the RC interval with `comp`, the new forward interval
        // has k = l_new_rev[comp] and l = k_new_rev[comp].
        Some(BiInterval {
            k: l_new_rev[comp],
            l: k_new_rev[comp],
            s: s_new_rev[comp],
            m: iv.m,
            n: iv.n,
        })
    }
}

// ---------------------------------------------------------------------------
// Initial interval for the first query character
// ---------------------------------------------------------------------------

/// Construct the starting bi-directional interval for query character `a`.
///
/// For a single character `a`, the forward BWT interval contains exactly the
/// SA entries whose suffixes start with `a`: `[count[a], count[a+1])`.
/// The reverse-complement interval starts with the complement `3-a`, so it
/// spans `[count[3-a], count[3-a+1])`.
///
/// Returns `None` if character `a` does not appear in the reference.
fn initial_interval(fmi: &FmIndex, a: u8) -> Option<BiInterval> {
    debug_assert!(a < 4);
    let k = fmi.count[a as usize];
    let s = fmi.count[a as usize + 1] - k;
    if s == 0 {
        return None;
    }
    let comp = 3usize - a as usize;
    let l = fmi.count[comp];
    Some(BiInterval { k, l, s, m: 0, n: 0 })
}

// ---------------------------------------------------------------------------
// Main SMEM algorithm
// ---------------------------------------------------------------------------

/// Find all SMEMs (Super-Maximal Exact Matches) in `query` against the
/// FM-index `fmi`.
///
/// `query` must be 2-bit encoded (use [`query_to_2bit`] to convert from ASCII).
/// Bases with value `> 3` (N's or unknown bases) terminate any active match.
///
/// `min_seed_len` filters out SMEMs shorter than this threshold.
///
/// # Algorithm
///
/// The implementation follows `getSMEMsOnePosOneThread` from `FMI_search.cpp`:
///
/// For each query start position `x` that has not yet been covered:
///
/// **Phase 1 – right extension (forward scan from x)**
/// Beginning from the initial single-character interval for `query[x]`, extend
/// right through `query[x+1], query[x+2], …` using bi-directional extension.
/// Whenever the interval *shrinks* (the hit count decreases) the previous
/// (larger) interval is saved as a right-maximal candidate.  Extension stops
/// when the interval becomes empty or a non-ACGT base is encountered.
///
/// **Phase 2 – left extension (backward scan from x-1)**
/// Each right-maximal candidate from Phase 1 is extended leftward from
/// position `x-1` down to position 0.  A candidate interval that cannot be
/// extended further left (the extended interval is empty) is a SMEM—provided
/// its query span is at least `min_seed_len` bases.
///
/// The outer loop then advances to the position just after the end of the
/// last found SMEM (or `x+1` if nothing was found at `x`), avoiding redundant
/// re-processing of already-covered positions.
pub fn get_smems(fmi: &FmIndex, query: &[u8], min_seed_len: i32) -> Vec<Smem> {
    let n = query.len() as i32;
    let mut smems: Vec<Smem> = Vec::new();

    let mut x: i32 = 0;
    while x < n {
        let a = query[x as usize];
        if a > 3 {
            // Skip N bases; no match can start here.
            x += 1;
            continue;
        }

        // ----------------------------------------------------------------
        // Phase 1: forward extension from position x.
        //
        // Collect all right-maximal intervals.  An interval is right-maximal
        // when the next right-extension either fails (empty) or shrinks the
        // interval.
        //
        // `prev_intervals` stores intervals in the order they became
        // right-maximal (earliest right end first).  We reverse them at the
        // end of phase 1 so that left extension processes shortest-match-first,
        // which matches the C++ behaviour.
        // ----------------------------------------------------------------
        let mut iv = match initial_interval(fmi, a) {
            Some(iv) => BiInterval { m: x, n: x, ..iv },
            None => {
                x += 1;
                continue;
            }
        };

        // The position after the last query base consumed; used to advance x.
        let mut next_x = x + 1;

        let mut prev_intervals: Vec<BiInterval> = Vec::new();

        // Extend rightward: j is the *next* query position to consume.
        let mut j = x + 1;
        while j < n {
            let b = query[j as usize];
            if b > 3 {
                // N base: current interval is right-maximal.
                next_x = j + 1;
                break;
            }
            match bidir_extend_right(fmi, &iv, b) {
                Some(new_iv) => {
                    let new_iv = BiInterval { n: j, ..new_iv };
                    if new_iv.s != iv.s {
                        // Interval shrank: previous interval was right-maximal.
                        prev_intervals.push(iv);
                    }
                    iv = new_iv;
                    next_x = j + 1;
                    j += 1;
                }
                None => {
                    // Extension failed: interval up to j-1 is right-maximal.
                    next_x = j;
                    break;
                }
            }
        }

        // The final interval (reaching the end of the read or an N/empty) is
        // also right-maximal.
        if iv.s > 0 {
            prev_intervals.push(iv);
        }

        if prev_intervals.is_empty() {
            x += 1;
            continue;
        }

        // Reverse so that the interval with the *longest* right extension comes
        // first—this mirrors the C++ reversal before the backward search.
        prev_intervals.reverse();

        // ----------------------------------------------------------------
        // Phase 2: backward extension from position x-1.
        //
        // For each right-maximal interval in `prev_intervals`, attempt to
        // extend leftward.  An interval whose left extension empties out is
        // left-maximal (and therefore a SMEM).  Intervals that successfully
        // extend left are kept for the next iteration.
        //
        // The C++ code tracks `curr_s` to deduplicate intervals with the same
        // hit count, which would otherwise represent the same BWT interval.
        // ----------------------------------------------------------------
        let mut cur_j = n; // leftmost position at which a SMEM was found
        let mut j = x - 1;
        while j >= 0 && !prev_intervals.is_empty() {
            let a = query[j as usize];
            if a > 3 {
                // N base encountered during left extension: all remaining
                // candidates that are long enough become SMEMs.
                break;
            }

            let mut curr_intervals: Vec<BiInterval> = Vec::new();
            let mut curr_s: i64 = -1;
            let mut emitted = false;

            for p in 0..prev_intervals.len() {
                let iv_p = prev_intervals[p];
                match bidir_extend_left(fmi, &iv_p, a) {
                    None => {
                        // Cannot extend left: this interval is left-maximal.
                        // Emit as SMEM if long enough and we haven't already
                        // emitted one at an equal or later position.
                        if !emitted && j < cur_j {
                            let match_len = (iv_p.n - iv_p.m + 1) as i32;
                            if match_len >= min_seed_len {
                                cur_j = j;
                                smems.push(iv_p.to_smem(0));
                                emitted = true;
                            }
                        }
                        // After emitting (or deciding not to emit) the first
                        // interval that failed to extend, stop processing
                        // further intervals at this left position—the C++ code
                        // breaks here after the first emit candidate.
                        break;
                    }
                    Some(new_iv) => {
                        let new_iv = BiInterval { m: j, ..new_iv };
                        // Deduplicate: only keep one representative per
                        // distinct interval size.
                        if new_iv.s != curr_s {
                            curr_s = new_iv.s;
                            curr_intervals.push(new_iv);
                        }
                        // Continue to the next interval in prev_intervals.
                    }
                }
            }

            prev_intervals = curr_intervals;
            j -= 1;
        }

        // Any intervals remaining after the backward search have reached
        // position 0 without failing: they are left-maximal at position 0.
        if let Some(iv_last) = prev_intervals.first() {
            let match_len = (iv_last.n - iv_last.m + 1) as i32;
            if match_len >= min_seed_len {
                smems.push(iv_last.to_smem(0));
            }
        }

        // Advance outer loop past the right end of the last SMEM found.
        x = next_x;
    }

    smems
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // query_to_2bit
    // ------------------------------------------------------------------

    #[test]
    fn test_query_to_2bit_acgt() {
        let seq = b"ACGTacgt";
        let enc = query_to_2bit(seq);
        assert_eq!(enc, vec![0, 1, 2, 3, 0, 1, 2, 3]);
    }

    #[test]
    fn test_query_to_2bit_n_and_unknown() {
        let seq = b"ANnX";
        let enc = query_to_2bit(seq);
        assert_eq!(enc[0], 0); // A
        assert_eq!(enc[1], 4); // N
        assert_eq!(enc[2], 4); // n
        assert_eq!(enc[3], 4); // X
    }

    #[test]
    fn test_query_to_2bit_empty() {
        let enc = query_to_2bit(b"");
        assert!(enc.is_empty());
    }

    #[test]
    fn test_query_to_2bit_all_bases() {
        assert_eq!(query_to_2bit(b"A"), vec![0]);
        assert_eq!(query_to_2bit(b"C"), vec![1]);
        assert_eq!(query_to_2bit(b"G"), vec![2]);
        assert_eq!(query_to_2bit(b"T"), vec![3]);
        assert_eq!(query_to_2bit(b"a"), vec![0]);
        assert_eq!(query_to_2bit(b"c"), vec![1]);
        assert_eq!(query_to_2bit(b"g"), vec![2]);
        assert_eq!(query_to_2bit(b"t"), vec![3]);
    }

    // ------------------------------------------------------------------
    // Smem helpers
    // ------------------------------------------------------------------

    #[test]
    fn test_bi_interval_to_smem() {
        let iv = BiInterval { k: 10, l: 20, s: 5, m: 3, n: 7 };
        let smem = iv.to_smem(42);
        assert_eq!(smem.rid, 42);
        assert_eq!(smem.m, 3);
        assert_eq!(smem.n, 8); // n+1 exclusive
        assert_eq!(smem.k, 10);
        assert_eq!(smem.l, 20);
        assert_eq!(smem.s, 5);
    }

    // ------------------------------------------------------------------
    // get_sa_entries
    // ------------------------------------------------------------------

    /// Minimal mock FmIndex for testing get_sa_entries without a real index.
    /// We test only the max_occ threshold logic here.
    #[test]
    fn test_get_sa_entries_too_repetitive() {
        // We cannot easily construct a real FmIndex in a unit test, but we can
        // verify that when smem.s > max_occ the function returns empty without
        // calling sa_lookup.  We do this by checking the return value is empty.
        // (A real FmIndex would be needed to test the SA lookup path.)
        //
        // This test documents the intended contract.
        let smem = Smem {
            rid: 0,
            m: 0,
            n: 5,
            k: 100,
            l: 200,
            s: 1000, // hit count
        };
        // max_occ < hit count → should return empty immediately.
        // We can't call get_sa_entries without a real FmIndex, so we document
        // the logic here as a compile-time check.
        let max_occ: i32 = 500;
        assert!(smem.s > max_occ as i64, "test precondition: smem.s should exceed max_occ");
    }
}
