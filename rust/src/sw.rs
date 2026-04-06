/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of banded Smith-Waterman alignment from bandedSWA.cpp / ksw.cpp.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

//! Banded Smith-Waterman local alignment.

// ---------------------------------------------------------------------------
// CIGAR constants and helpers (BAM encoding)
// ---------------------------------------------------------------------------

pub const CIGAR_M: u32 = 0; // match/mismatch
pub const CIGAR_I: u32 = 1; // insertion into reference
pub const CIGAR_D: u32 = 2; // deletion from reference
pub const CIGAR_S: u32 = 4; // soft clip

#[inline] pub fn cigar_op(c: u32)  -> u32 { c & 0xf }
#[inline] pub fn cigar_len(c: u32) -> u32 { c >> 4 }
#[inline] pub fn make_cigar(op: u32, len: u32) -> u32 { (len << 4) | (op & 0xf) }

// ---------------------------------------------------------------------------
// Result type
// ---------------------------------------------------------------------------

/// Result of a Smith-Waterman alignment.
#[derive(Debug, Clone, Default)]
pub struct SwResult {
    /// Best alignment score.
    pub score: i32,
    /// Second-best score (sub-optimal).
    pub score2: i32,
    /// Query begin (0-based, inclusive).
    pub qb: i32,
    /// Query end (0-based, exclusive).
    pub qe: i32,
    /// Reference begin (relative to the supplied reference slice).
    pub rb: i32,
    /// Reference end (relative to the supplied reference slice).
    pub re: i32,
    /// CIGAR ops in BAM encoding.
    pub cigar: Vec<u32>,
    /// Number of CIGAR ops.
    pub n_cigar: i32,
}

// ---------------------------------------------------------------------------
// Scoring helpers
// ---------------------------------------------------------------------------

/// Look up the score for aligning query base `q` vs reference base `r` using
/// a 5×5 scoring matrix (ACGTN, row-major: `mat[q*5 + r]`).
#[inline]
fn score_pair(mat: &[i8], m: usize, q: u8, r: u8) -> i32 {
    let qi = if (q as usize) < m { q as usize } else { m - 1 };
    let ri = if (r as usize) < m { r as usize } else { m - 1 };
    mat[qi * m + ri] as i32
}

// ---------------------------------------------------------------------------
// ksw_extend: banded extension without traceback
// ---------------------------------------------------------------------------

/// Banded extension from a seed.  Computes the best score and end positions
/// without producing a CIGAR string.
///
/// The band covers columns `[i - w, i + w]` for query row `i`.
/// Returns a `SwResult` with `score`, `qb`, `qe`, `rb`, `re` set but
/// `cigar` empty (use [`banded_sw`] if you need the CIGAR).
pub fn ksw_extend(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    zdrop: i32,
    _end_bonus: i32,
    h0: i32,
) -> SwResult {
    let qlen = qlen as usize;
    let tlen = tlen as usize;
    let m = m as usize;
    let w = w as usize;

    if qlen == 0 || tlen == 0 {
        return SwResult { qe: qlen as i32, re: tlen as i32, ..Default::default() };
    }

    // H[j] = best score ending at query pos i, ref pos j
    // E[j] = best score of gap in query ending at ref pos j (insertion into ref)
    let mut h = vec![i32::MIN / 2; tlen + 1];
    let mut e = vec![i32::MIN / 2; tlen + 1];

    // Initialise first row using h0
    h[0] = h0;
    for j in 1..=tlen {
        h[j] = h[j - 1] - if j == 1 { o_ins + e_ins } else { e_ins };
        if h[j] < 0 { h[j] = 0; }
    }

    let mut best = h0;
    let mut best_i = 0i32;
    let mut best_j = 0i32;
    let mut qb = 0i32;
    let mut rb = 0i32;

    for i in 0..qlen {
        let q_base = query[i];
        let mut f = i32::MIN / 2;
        let mut h_prev = if i > 0 { 0 } else { h0 };
        let j_lo = if i >= w { i - w } else { 0 };
        let j_hi = (i + w + 1).min(tlen);

        let mut prev_h = h[j_lo];

        for j in j_lo..j_hi {
            let t_base = target[j];
            let diag_score = score_pair(mat, m, q_base, t_base);

            // H(i,j) = max(H(i-1,j-1)+score, E(i,j), F(i,j))
            let h_diag = if j > 0 { h[j - 1] } else { h_prev };
            let h_cur = (h_diag + diag_score)
                .max(e[j])
                .max(f);
            let h_cur = h_cur.max(0);

            // Update E (gap in query = insertion into ref)
            e[j] = (h_cur - o_ins - e_ins).max(e[j] - e_ins);
            // Update F (gap in ref = deletion from ref)
            f = (h_cur - o_del - e_del).max(f - e_del);

            prev_h = h[j];
            h[j] = h_cur;

            if h_cur > best {
                best = h_cur;
                best_i = i as i32;
                best_j = j as i32;
            }
        }
        h_prev = prev_h;

        // Z-drop: if score has dropped too much from the best, terminate
        if zdrop > 0 && best - h[j_hi.saturating_sub(1)] > zdrop {
            break;
        }
        let _ = h_prev;
    }

    SwResult {
        score: best,
        qb,
        qe: best_i + 1,
        rb,
        re: best_j + 1,
        ..Default::default()
    }
}

// ---------------------------------------------------------------------------
// banded_sw: full banded DP with traceback
// ---------------------------------------------------------------------------

/// Direction flags for traceback.
const DIR_DIAG: u8 = 0;
const DIR_INS:  u8 = 1; // gap in ref (F matrix)
const DIR_DEL:  u8 = 2; // gap in query (E matrix)

/// Full banded Smith-Waterman alignment with CIGAR traceback.
///
/// `query` and `reference` are 2-bit encoded (A=0,C=1,G=2,T=3,N=4).
/// `score_init` is the starting score (usually 0).
pub fn banded_sw(
    query: &[u8],
    reference: &[u8],
    qlen: i32,
    rlen: i32,
    score_init: i32,
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
    mat: &[i8],
    m: i32,
) -> SwResult {
    let qlen = qlen as usize;
    let rlen = rlen as usize;
    let w = w as usize;
    let m = m as usize;

    if qlen == 0 || rlen == 0 {
        return SwResult { qe: qlen as i32, re: rlen as i32, ..Default::default() };
    }

    // Allocate H, E, F and direction matrices
    // We use a flat array indexed [i*(rlen+1) + j]
    let rows = qlen + 1;
    let cols = rlen + 1;
    let mut h_mat = vec![0i32; rows * cols];
    let mut e_mat = vec![i32::MIN / 2; rows * cols];
    let mut f_mat = vec![i32::MIN / 2; rows * cols];
    let mut dir   = vec![DIR_DIAG; rows * cols];

    // Initialize
    h_mat[0] = score_init;

    // Fill first row: gaps in query (E matrix)
    for j in 1..cols {
        let prev = h_mat[j - 1];
        let gap_cost = if j == 1 { o_del + e_del } else { e_del };
        h_mat[j] = (prev - gap_cost).max(0);
        dir[j] = DIR_DEL;
    }
    // Fill first column: gaps in ref (F matrix)
    for i in 1..rows {
        let prev = h_mat[(i - 1) * cols];
        let gap_cost = if i == 1 { o_ins + e_ins } else { e_ins };
        h_mat[i * cols] = (prev - gap_cost).max(0);
        dir[i * cols] = DIR_INS;
    }

    let mut best_score = score_init;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..rows {
        let q_base = query[i - 1];
        let j_lo = if i > w + 1 { i - w - 1 } else { 1 };
        let j_hi = (i + w).min(rlen) + 1;

        for j in j_lo..j_hi {
            let r_base = reference[j - 1];
            let idx = i * cols + j;

            let sc = score_pair(mat, m, q_base, r_base);
            let h_diag = h_mat[(i - 1) * cols + (j - 1)];

            // E: gap in query (deletion from ref) – extend left
            let e_open  = h_mat[i * cols + (j - 1)] - o_del - e_del;
            let e_ext   = e_mat[i * cols + (j - 1)] - e_del;
            let e_val   = e_open.max(e_ext);
            e_mat[idx]  = e_val;

            // F: gap in ref (insertion into ref) – extend up
            let f_open  = h_mat[(i - 1) * cols + j] - o_ins - e_ins;
            let f_ext   = f_mat[(i - 1) * cols + j] - e_ins;
            let f_val   = f_open.max(f_ext);
            f_mat[idx]  = f_val;

            // H = max(diag + score, E, F, 0)
            let mut h_val = h_diag + sc;
            let mut d = DIR_DIAG;
            if e_val > h_val { h_val = e_val; d = DIR_DEL; }
            if f_val > h_val { h_val = f_val; d = DIR_INS; }
            if h_val < 0 { h_val = 0; d = DIR_DIAG; }
            h_mat[idx] = h_val;
            dir[idx] = d;

            if h_val > best_score {
                best_score = h_val;
                best_i = i;
                best_j = j;
            }
        }
    }

    if best_score <= 0 {
        return SwResult::default();
    }

    // Traceback from (best_i, best_j)
    let mut cigar_ops: Vec<(u32, u32)> = Vec::new(); // (op, count)
    let (mut ti, mut tj) = (best_i, best_j);

    while ti > 0 && tj > 0 && h_mat[ti * cols + tj] > 0 {
        let d = dir[ti * cols + tj];
        let op = match d {
            DIR_DIAG => { ti -= 1; tj -= 1; CIGAR_M }
            DIR_INS  => { ti -= 1; CIGAR_I }
            _        => { tj -= 1; CIGAR_D }
        };
        match cigar_ops.last_mut() {
            Some(last) if last.0 == op => last.1 += 1,
            _ => cigar_ops.push((op, 1)),
        }
    }

    cigar_ops.reverse();
    let cigar: Vec<u32> = cigar_ops.iter().map(|&(op, len)| make_cigar(op, len)).collect();
    let n_cigar = cigar.len() as i32;

    SwResult {
        score: best_score,
        score2: 0,
        qb: ti as i32,
        qe: best_i as i32,
        rb: tj as i32,
        re: best_j as i32,
        cigar,
        n_cigar,
    }
}

// ---------------------------------------------------------------------------
// ksw_global: global (Needleman-Wunsch) banded alignment
// ---------------------------------------------------------------------------

/// Banded global alignment (Needleman-Wunsch).  Used for the final alignment
/// extension where start and end are fixed.
pub fn ksw_global(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    w: i32,
) -> SwResult {
    let qlen = qlen as usize;
    let tlen = tlen as usize;
    let w = w as usize;
    let m = m as usize;

    if qlen == 0 || tlen == 0 {
        return SwResult::default();
    }

    let rows = qlen + 1;
    let cols = tlen + 1;
    let mut h  = vec![i32::MIN / 2; rows * cols];
    let mut e  = vec![i32::MIN / 2; rows * cols];
    let mut f  = vec![i32::MIN / 2; rows * cols];
    let mut dir = vec![DIR_DIAG; rows * cols];

    h[0] = 0;
    for j in 1..cols {
        h[j] = -(o_del + e_del * j as i32);
        dir[j] = DIR_DEL;
    }
    for i in 1..rows {
        h[i * cols] = -(o_ins + e_ins * i as i32);
        dir[i * cols] = DIR_INS;
    }

    for i in 1..rows {
        let q_base = query[i - 1];
        let j_lo = if i > w + 1 { i - w - 1 } else { 1 };
        let j_hi = (i + w).min(tlen) + 1;

        for j in j_lo..j_hi {
            let r_base = target[j - 1];
            let idx = i * cols + j;
            let sc = score_pair(mat, m, q_base, r_base);

            let e_val = (h[i * cols + j - 1] - o_del - e_del).max(e[i * cols + j - 1] - e_del);
            e[idx] = e_val;

            let f_val = (h[(i - 1) * cols + j] - o_ins - e_ins).max(f[(i - 1) * cols + j] - e_ins);
            f[idx] = f_val;

            let mut h_val = h[(i - 1) * cols + j - 1] + sc;
            let mut d = DIR_DIAG;
            if e_val > h_val { h_val = e_val; d = DIR_DEL; }
            if f_val > h_val { h_val = f_val; d = DIR_INS; }
            h[idx] = h_val;
            dir[idx] = d;
        }
    }

    let best_score = h[rows * cols - 1];

    // Traceback
    let mut cigar_ops: Vec<(u32, u32)> = Vec::new();
    let (mut ti, mut tj) = (qlen, tlen);
    while ti > 0 || tj > 0 {
        if ti == 0 { cigar_ops.push((CIGAR_D, tj as u32)); break; }
        if tj == 0 { cigar_ops.push((CIGAR_I, ti as u32)); break; }
        let d = dir[ti * cols + tj];
        let op = match d {
            DIR_DIAG => { ti -= 1; tj -= 1; CIGAR_M }
            DIR_INS  => { ti -= 1; CIGAR_I }
            _        => { tj -= 1; CIGAR_D }
        };
        match cigar_ops.last_mut() {
            Some(last) if last.0 == op => last.1 += 1,
            _ => cigar_ops.push((op, 1)),
        }
    }
    cigar_ops.reverse();
    let cigar: Vec<u32> = cigar_ops.iter().map(|&(op, len)| make_cigar(op, len)).collect();
    let n_cigar = cigar.len() as i32;

    SwResult { score: best_score, qb: 0, qe: qlen as i32, rb: 0, re: tlen as i32, cigar, n_cigar, ..Default::default() }
}

// ---------------------------------------------------------------------------
// ksw_align: local SW (Smith-Waterman) without banding
// ---------------------------------------------------------------------------

/// Standard Smith-Waterman local alignment.
/// Returns score, end positions, and CIGAR (via traceback).
pub fn ksw_align(
    qlen: i32,
    query: &[u8],
    tlen: i32,
    target: &[u8],
    m: i32,
    mat: &[i8],
    o_del: i32,
    e_del: i32,
    o_ins: i32,
    e_ins: i32,
    _xtra: i32,
) -> SwResult {
    // Use a large band width so all cells are computed
    let w = qlen.max(tlen) + 1;
    banded_sw(query, target, qlen, tlen, 0, o_del, e_del, o_ins, e_ins, w, mat, m)
}

// ---------------------------------------------------------------------------
// Edit-distance from CIGAR
// ---------------------------------------------------------------------------

/// Compute edit distance (NM tag) from a CIGAR and the aligned sequences.
pub fn compute_nm_cigar(
    query: &[u8],
    target: &[u8],
    cigar: &[u32],
    qb: i32,
    rb: i32,
) -> i32 {
    let mut nm = 0i32;
    let mut qi = qb as usize;
    let mut ri = rb as usize;
    for &c in cigar {
        let op  = cigar_op(c);
        let len = cigar_len(c) as usize;
        match op {
            o if o == CIGAR_M => {
                for k in 0..len {
                    let qb = if qi + k < query.len()  { query[qi + k]  } else { 4 };
                    let rb = if ri + k < target.len() { target[ri + k] } else { 4 };
                    if qb != rb { nm += 1; }
                }
                qi += len;
                ri += len;
            }
            o if o == CIGAR_I => { nm += len as i32; qi += len; }
            o if o == CIGAR_D => { nm += len as i32; ri += len; }
            o if o == CIGAR_S => { qi += len; }
            _ => {}
        }
    }
    nm
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_mat() -> [i8; 25] {
        let mut m = [0i8; 25];
        // Match=2, mismatch=-3
        for i in 0..4 {
            for j in 0..4 {
                m[i * 5 + j] = if i == j { 2 } else { -3 };
            }
            m[i * 5 + 4] = -1;
        }
        for j in 0..5 { m[4 * 5 + j] = -1; }
        m
    }

    #[test]
    fn cigar_ops_basic() {
        let c = make_cigar(CIGAR_M, 10);
        assert_eq!(cigar_op(c), CIGAR_M);
        assert_eq!(cigar_len(c), 10);
    }

    #[test]
    fn banded_sw_identical() {
        let mat = simple_mat();
        // query = target = ACGT (2-bit)
        let seq: Vec<u8> = vec![0, 1, 2, 3];
        let res = banded_sw(&seq, &seq, 4, 4, 0, 6, 1, 6, 1, 4, &mat, 5);
        // Perfect match: 4 * 2 = 8
        assert_eq!(res.score, 8);
        assert!(!res.cigar.is_empty());
        let op = cigar_op(res.cigar[0]);
        assert_eq!(op, CIGAR_M);
    }

    #[test]
    fn banded_sw_empty() {
        let mat = simple_mat();
        let res = banded_sw(&[], &[], 0, 0, 0, 6, 1, 6, 1, 4, &mat, 5);
        assert_eq!(res.score, 0);
    }

    #[test]
    fn ksw_global_identical() {
        let mat = simple_mat();
        let seq: Vec<u8> = vec![0, 1, 2, 3];
        let res = ksw_global(4, &seq, 4, &seq, 5, &mat, 6, 1, 6, 1, 4);
        assert_eq!(res.score, 8); // 4 matches × 2
    }

    #[test]
    fn compute_nm_all_matches() {
        let q: Vec<u8> = vec![0, 1, 2, 3];
        let r: Vec<u8> = vec![0, 1, 2, 3];
        let cigar = vec![make_cigar(CIGAR_M, 4)];
        assert_eq!(compute_nm_cigar(&q, &r, &cigar, 0, 0), 0);
    }

    #[test]
    fn compute_nm_one_mismatch() {
        let q: Vec<u8> = vec![0, 1, 2, 3];
        let r: Vec<u8> = vec![0, 1, 2, 0]; // last base differs
        let cigar = vec![make_cigar(CIGAR_M, 4)];
        assert_eq!(compute_nm_cigar(&q, &r, &cigar, 0, 0), 1);
    }
}
