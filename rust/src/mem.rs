/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of the core MEM alignment algorithm from bwamem.cpp.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

//! Core MEM (Maximal Exact Match) algorithm: seeding, chaining, extension.

use crate::bntseq::{BntSeq, pac2nt, restore_pac_rev};
use crate::bwt::FmIndex;
use crate::fmi_search::{get_sa_entries, get_smems, query_to_2bit};
use crate::sw::{banded_sw, ksw_extend, make_cigar, CIGAR_D, CIGAR_I, CIGAR_M, CIGAR_S};
use crate::types::{BSeq, MemAlnReg, MemChain, MemOpt, MemSeed};

// ---------------------------------------------------------------------------
// Seed collection
// ---------------------------------------------------------------------------

/// Collect MEM seeds for a query sequence.
///
/// 1. Run SMEM finding via the FM-index.
/// 2. Convert each SMEM to one or more `MemSeed` records via SA lookup.
/// 3. Filter seeds by valid reference coordinates.
pub fn mem_collect_seeds(
    opt: &MemOpt,
    fmi: &FmIndex,
    bns: &BntSeq,
    query_2bit: &[u8],
) -> Vec<MemSeed> {
    let smems = get_smems(fmi, query_2bit, opt.min_seed_len);
    let mut seeds: Vec<MemSeed> = Vec::new();

    for smem in &smems {
        let match_len = (smem.n - smem.m) as i32;
        let sa_entries = get_sa_entries(fmi, smem, opt.max_occ);

        for ref_pos in sa_entries {
            // The SA gives a position in the doubled reference [0, 2*l_pac).
            // Filter out positions that straddle the forward/reverse boundary.
            let rbeg = ref_pos;
            let rend = rbeg + match_len as i64;
            if rbeg < 0 { continue; }
            // Skip if the seed crosses the l_pac boundary
            if rbeg < bns.l_pac && rend > bns.l_pac { continue; }
            if rbeg >= 2 * bns.l_pac { continue; }

            seeds.push(MemSeed {
                rbeg,
                qbeg: smem.m as i32,
                len: match_len,
                score: opt.a * match_len,
                done: false,
            });
        }
    }

    // Sort by reference position, then query position
    seeds.sort_by(|a, b| a.rbeg.cmp(&b.rbeg).then(a.qbeg.cmp(&b.qbeg)));
    seeds
}

// ---------------------------------------------------------------------------
// Chaining
// ---------------------------------------------------------------------------

/// Chain seeds into colinear groups.
///
/// Seeds are chained if they are colinear (same strand, increasing positions)
/// and the gap between them does not exceed `opt.max_chain_gap`.
pub fn mem_chain_seeds(opt: &MemOpt, bns: &BntSeq, seeds: Vec<MemSeed>) -> Vec<MemChain> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let mut chains: Vec<MemChain> = Vec::new();

    'outer: for seed in seeds {
        // Try to extend an existing chain
        let mut best_chain: Option<usize> = None;
        let mut best_gap = i64::MAX;

        for (ci, chain) in chains.iter().enumerate() {
            if chain.kept == 0 { continue; }
            let last = chain.seeds.last().unwrap();

            // Must be on the same strand (both < l_pac or both >= l_pac)
            let same_strand = (seed.rbeg < bns.l_pac) == (last.rbeg < bns.l_pac);
            if !same_strand { continue; }

            // Must be colinear
            if seed.rbeg <= last.rbeg { continue; }
            if seed.qbeg <= last.qbeg { continue; }

            let rdiff = seed.rbeg - last.rbeg - last.len as i64;
            let qdiff = seed.qbeg - last.qbeg - last.len;

            // Gap too large
            if rdiff > opt.max_chain_gap as i64 { continue; }
            if qdiff > opt.max_chain_gap { continue; }
            if rdiff < -(last.len as i64) { continue; }
            if qdiff < -(last.len)       { continue; }

            // Prefer the chain with the smallest gap
            let gap = (rdiff - qdiff as i64).abs();
            if gap < best_gap {
                best_gap = gap;
                best_chain = Some(ci);
            }
        }

        if let Some(ci) = best_chain {
            let w = chains[ci].seeds.last().map(|s| s.len).unwrap_or(0);
            let seed_len = seed.len;
            chains[ci].seeds.push(seed);
            chains[ci].weight += seed_len;
            chains[ci].w = w.max(seed_len);
        } else {
            // Start a new chain
            let rid = bns.pos2rid(if seed.rbeg < bns.l_pac { seed.rbeg } else { 2 * bns.l_pac - 1 - seed.rbeg });
            let w = seed.len;
            let weight = seed.len;
            chains.push(MemChain {
                seqid: rid,
                seeds: vec![seed],
                w,
                kept: 1,
                is_alt: false,
                frac_rep: 0.0,
                weight,
            });
        }
    }

    chains
}

// ---------------------------------------------------------------------------
// Chain filtering
// ---------------------------------------------------------------------------

/// Filter chains:
/// 1. Compute weight for each chain.
/// 2. Remove chains with weight below min_chain_weight.
/// 3. Mark redundant chains (highly overlapping with a better chain).
/// 4. Limit to max_chain_extend chains.
pub fn mem_chain_flt(opt: &MemOpt, mut chains: Vec<MemChain>) -> Vec<MemChain> {
    // Sort by weight descending
    chains.sort_by(|a, b| b.weight.cmp(&a.weight));

    // Recompute weights and filter by min
    chains.retain(|c| c.weight >= opt.min_chain_weight);

    // Mark redundant chains
    let n = chains.len();
    for i in 1..n {
        let qi_beg = chains[i].seeds.first().map(|s| s.qbeg).unwrap_or(0);
        let qi_end = chains[i].seeds.last().map(|s| s.qbeg + s.len).unwrap_or(0);
        let qi_len = (qi_end - qi_beg).max(1);

        for j in 0..i {
            if chains[j].kept == 0 { continue; }
            let qj_beg = chains[j].seeds.first().map(|s| s.qbeg).unwrap_or(0);
            let qj_end = chains[j].seeds.last().map(|s| s.qbeg + s.len).unwrap_or(0);

            // Compute overlap
            let ovlp = (qi_end.min(qj_end) - qi_beg.max(qj_beg)).max(0);
            if ovlp as f32 / qi_len as f32 > opt.mask_level {
                // Chain i is largely covered by chain j
                chains[i].kept = 0;
                break;
            }
        }
    }

    // Retain only kept chains, up to max_chain_extend
    chains.retain(|c| c.kept != 0);
    chains.truncate(opt.max_chain_extend as usize);
    chains
}

// ---------------------------------------------------------------------------
// Reference extraction
// ---------------------------------------------------------------------------

/// Extract reference bases (2-bit encoded) for interval [beg, end) from the
/// doubled packed reference (fwd+rev).
fn get_ref_bases(pac: &[u8], l_pac: i64, beg: i64, end: i64) -> Vec<u8> {
    let len = (end - beg).max(0) as usize;
    let mut ref_seq = Vec::with_capacity(len);
    let total = 2 * l_pac;
    for pos in beg..end {
        let p = if pos >= 0 && pos < total { pos } else { continue };
        ref_seq.push(pac2nt(pac, p));
    }
    ref_seq
}

// ---------------------------------------------------------------------------
// Chain to alignment regions
// ---------------------------------------------------------------------------

/// Extend one chain into one or more alignment regions using banded SW.
pub fn mem_chain2aln(
    opt: &MemOpt,
    bns: &BntSeq,
    pac: &[u8],
    query: &[u8],
    query_2bit: &[u8],
    chain: &MemChain,
) -> Vec<MemAlnReg> {
    if chain.seeds.is_empty() {
        return Vec::new();
    }

    let qlen = query.len() as i32;
    let mut regs: Vec<MemAlnReg> = Vec::new();

    for seed in &chain.seeds {
        if seed.done { continue; }

        // Determine the band width
        let w = opt.w.max(seed.len / 2);

        // Reference region to align: extend by band width around the seed
        let ref_beg = (seed.rbeg - w as i64).max(0);
        let ref_end = (seed.rbeg + seed.len as i64 + w as i64).min(2 * bns.l_pac);

        // Query region: with some padding
        let q_beg = (seed.qbeg - w).max(0) as usize;
        let q_end = ((seed.qbeg + seed.len + w) as usize).min(query.len());

        if q_beg >= q_end || ref_beg >= ref_end { continue; }

        let ref_seq = get_ref_bases(pac, bns.l_pac, ref_beg, ref_end);
        let qslice = &query_2bit[q_beg..q_end];

        let sw_res = banded_sw(
            qslice,
            &ref_seq,
            (q_end - q_beg) as i32,
            ref_seq.len() as i32,
            0,
            opt.o_del,
            opt.e_del,
            opt.o_ins,
            opt.e_ins,
            w,
            &opt.mat,
            5,
        );

        if sw_res.score < opt.t { continue; }

        let rb = ref_beg + sw_res.rb as i64;
        let re = ref_beg + sw_res.re as i64;
        let qb = q_beg as i32 + sw_res.qb;
        let qe = q_beg as i32 + sw_res.qe;

        let rid = bns.pos2rid(if rb < bns.l_pac { rb } else { 2 * bns.l_pac - 1 - rb });

        let mut reg = MemAlnReg {
            rb,
            re,
            qb,
            qe,
            rid,
            score: sw_res.score,
            truesc: sw_res.score,
            sub: 0,
            alt_sc: 0,
            csub: 0,
            sub_n: 0,
            w,
            seedcov: seed.len,
            secondary: -1,
            secondary_all: -1,
            seedlen0: seed.len,
            n_comp: 1,
            is_alt: false,
            frac_rep: chain.frac_rep,
            hash: 0,
            cigar: sw_res.cigar,
            n_cigar: sw_res.n_cigar,
            mapq: 0,
        };

        regs.push(reg);
    }

    regs
}

// ---------------------------------------------------------------------------
// Primary/secondary marking
// ---------------------------------------------------------------------------

/// Mark primary vs secondary alignments.
/// A region is secondary if it significantly overlaps a higher-scoring one.
pub fn mem_mark_primary(opt: &MemOpt, regs: &mut Vec<MemAlnReg>) {
    // Sort by score descending
    regs.sort_by(|a, b| b.score.cmp(&a.score));

    let n = regs.len();
    for i in 0..n {
        let (qi_beg, qi_end, ri_beg, ri_end) = {
            let r = &regs[i];
            (r.qb, r.qe, r.rb, r.re)
        };
        let qi_len = (qi_end - qi_beg).max(1);

        for j in 0..i {
            if regs[j].secondary >= 0 { continue; }
            let (qj_beg, qj_end) = (regs[j].qb, regs[j].qe);
            let (rj_beg, rj_end) = (regs[j].rb, regs[j].re);

            // Check query overlap
            let q_ovlp = (qi_end.min(qj_end) - qi_beg.max(qj_beg)).max(0) as f32;
            // Check reference overlap
            let r_ovlp = ((ri_end.min(rj_end) - ri_beg.max(rj_beg)).max(0)) as f32;
            let r_len = (ri_end - ri_beg).max(1) as f32;

            if q_ovlp / qi_len as f32 > opt.mask_level || r_ovlp / r_len > opt.mask_level {
                regs[i].secondary = j as i32;
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MAPQ estimation
// ---------------------------------------------------------------------------

/// Estimate mapping quality for a primary alignment.
pub fn mem_approx_mapq(opt: &MemOpt, reg: &MemAlnReg) -> u8 {
    if reg.secondary >= 0 {
        return 0;
    }
    if reg.sub <= 0 {
        return 60;
    }
    let score = reg.score as f64;
    let sub   = reg.sub   as f64;
    let a     = opt.a     as f64;
    // MAPQ formula from bwa-mem: 6.02 * (score-sub) / a + 0.499
    let raw = 6.02 * (score - sub) / a + 0.499;
    (raw as u8).min(60)
}

// ---------------------------------------------------------------------------
// Per-read alignment pipeline
// ---------------------------------------------------------------------------

/// Align one read end-to-end and return all alignment regions.
pub fn mem_align1(
    opt: &MemOpt,
    fmi: &FmIndex,
    bns: &BntSeq,
    pac: &[u8],
    seq: &BSeq,
) -> Vec<MemAlnReg> {
    if seq.seq.is_empty() {
        return Vec::new();
    }

    // Convert to 2-bit
    let query_2bit = query_to_2bit(&seq.seq);

    // 1. Collect seeds
    let seeds = mem_collect_seeds(opt, fmi, bns, &query_2bit);

    // 2. Chain seeds
    let chains = mem_chain_seeds(opt, bns, seeds);

    // 3. Filter chains
    let chains = mem_chain_flt(opt, chains);

    // 4. Extend chains to alignment regions
    let mut regs: Vec<MemAlnReg> = Vec::new();
    for chain in &chains {
        let mut chain_regs = mem_chain2aln(opt, bns, pac, &seq.seq, &query_2bit, chain);
        // Set MAPQ placeholder
        for r in &mut chain_regs {
            r.mapq = 0;
        }
        regs.extend(chain_regs);
    }

    // 5. Mark primary/secondary
    mem_mark_primary(opt, &mut regs);

    // 6. Compute MAPQ
    let n = regs.len();
    // Find best sub score for primary MAPQ
    let best_score = regs.iter().filter(|r| r.secondary < 0).map(|r| r.score).max().unwrap_or(0);
    let sub_score = regs.iter().filter(|r| r.secondary < 0)
        .map(|r| r.score).filter(|&s| s < best_score).max().unwrap_or(0);

    for reg in &mut regs {
        if reg.secondary < 0 {
            reg.sub = sub_score;
            reg.mapq = mem_approx_mapq(opt, reg);
        }
    }

    regs
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::MemOpt;

    #[test]
    fn chain_flt_min_weight() {
        let opt = MemOpt { min_chain_weight: 20, ..MemOpt::default() };
        let chains = vec![
            MemChain { weight: 10, kept: 1, ..Default::default() },
            MemChain { weight: 30, kept: 1, ..Default::default() },
        ];
        let flt = mem_chain_flt(&opt, chains);
        assert_eq!(flt.len(), 1);
        assert_eq!(flt[0].weight, 30);
    }

    #[test]
    fn mem_approx_mapq_no_sub() {
        let opt = MemOpt::default();
        let reg = MemAlnReg { score: 100, sub: 0, secondary: -1, ..Default::default() };
        assert_eq!(mem_approx_mapq(&opt, &reg), 60);
    }

    #[test]
    fn mem_approx_mapq_secondary() {
        let opt = MemOpt::default();
        let reg = MemAlnReg { score: 100, sub: 90, secondary: 0, ..Default::default() };
        assert_eq!(mem_approx_mapq(&opt, &reg), 0);
    }
}
