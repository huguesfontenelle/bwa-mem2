/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of core data structures from bwamem.h / bwa.h.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

// ---------------------------------------------------------------------------
// MEM_F_* flag constants  (from bwamem.h)
// ---------------------------------------------------------------------------

pub const MEM_F_PE: i32           = 0x2;
pub const MEM_F_NOPAIRING: i32    = 0x4;
pub const MEM_F_ALL: i32          = 0x8;
pub const MEM_F_NO_MULTI: i32     = 0x10;
pub const MEM_F_NO_RESCUE: i32    = 0x20;
pub const MEM_F_REF_HDR: i32      = 0x100;
pub const MEM_F_SOFTCLIP: i32     = 0x200;
pub const MEM_F_SMARTPE: i32      = 0x400;
pub const MEM_F_PRIMARY5: i32     = 0x800;
pub const MEM_F_KEEP_SUPP_MAPQ: i32 = 0x1000;
pub const MEM_F_XB: i32           = 0x2000;

pub const MEM_MAPQ_COEF: f64 = 30.0;
pub const MEM_MAPQ_MAX: i32  = 60;

// ---------------------------------------------------------------------------
// Scoring / alignment options  (mem_opt_t in C)
// ---------------------------------------------------------------------------

/// All tuning parameters for the BWA-MEM2 aligner, mirroring `mem_opt_t`.
#[derive(Debug, Clone)]
pub struct MemOpt {
    /// Match score (positive).
    pub a: i32,
    /// Mismatch penalty (positive value; applied as -b internally).
    pub b: i32,
    /// Gap-open penalty for deletions.
    pub o_del: i32,
    /// Gap-extension penalty for deletions.
    pub e_del: i32,
    /// Gap-open penalty for insertions.
    pub o_ins: i32,
    /// Gap-extension penalty for insertions.
    pub e_ins: i32,
    /// Phred-scaled penalty for unpaired reads in PE mode.
    pub pen_unpaired: i32,
    /// Soft/hard clipping penalty at the 5′ end.
    pub pen_clip5: i32,
    /// Soft/hard clipping penalty at the 3′ end.
    pub pen_clip3: i32,
    /// Band width for banded Smith-Waterman.
    pub w: i32,
    /// Z-dropoff threshold; terminate extension when score drops by this much.
    pub zdrop: i32,
    /// Maximum seed interval size considered for exact-match extension.
    pub max_mem_intv: i64,
    /// Output score threshold (hits below this are discarded).
    pub t: i32,
    /// Combination of `MEM_F_*` flags.
    pub flag: i32,
    /// Minimum seed length.
    pub min_seed_len: i32,
    /// Minimum chain weight (seed coverage) to keep a chain.
    pub min_chain_weight: i32,
    /// Maximum number of seeds to extend per chain.
    pub max_chain_extend: i32,
    /// Split a seed into two if the MEM length exceeds `min_seed_len * split_factor`.
    pub split_factor: f32,
    /// Split a seed if its occurrence count falls below this value.
    pub split_width: i32,
    /// Skip a seed if its occurrence count exceeds this value.
    pub max_occ: i32,
    /// Maximum gap between seeds in the same chain.
    pub max_chain_gap: i32,
    /// Number of worker threads.
    pub n_threads: i32,
    /// Number of bases processed per batch.
    pub chunk_size: i64,
    /// A hit is considered redundant if its overlap with a better hit exceeds
    /// `mask_level` times the shorter hit length.
    pub mask_level: f32,
    /// Drop a chain whose seed coverage is below `drop_ratio` times the best
    /// overlapping chain's coverage.
    pub drop_ratio: f32,
    /// When building the XA tag, ignore hits with score < `xa_drop_ratio * max_score`.
    pub xa_drop_ratio: f32,
    pub mask_level_redun: f32,
    pub mapq_coeff_len: f32,
    /// `ln(mapq_coeff_len)`, precomputed.
    pub mapq_coeff_fac: f32,
    /// Skip PE pairs with insert size larger than this when estimating insert-size dist.
    pub max_ins: i32,
    /// Maximum rounds of mate-SW rescue per end.
    pub max_matesw: i32,
    /// Maximum number of XA hits reported for primary-chromosome hits.
    pub max_xs_aux: i32,
    /// (Unused in current code; reserved for min insert size filtering.)
    pub min_ins: i32,
    /// 5×5 Smith-Waterman scoring matrix for {A,C,G,T,N}.
    /// Stored row-major: `mat[i*5 + j]` is the score for aligning base `i` to base `j`.
    pub mat: [i8; 25],
}

/// Fill a 5×5 scoring matrix from match score `a` and mismatch penalty `b`.
///
/// Matches `bwa_fill_scmat` in `bwa.cpp`.
pub fn fill_scmat(a: i32, b: i32, mat: &mut [i8; 25]) {
    for i in 0..4usize {
        for j in 0..4usize {
            mat[i * 5 + j] = if i == j { a as i8 } else { -(b as i8) };
        }
        mat[i * 5 + 4] = -1_i8; // N column
    }
    // N row
    for j in 0..5usize {
        mat[4 * 5 + j] = -1_i8;
    }
}

impl Default for MemOpt {
    /// Initialise with the same defaults as `mem_opt_init()` in `bwamem.cpp`.
    fn default() -> Self {
        let a: i32 = 1;
        let b: i32 = 4;
        let mapq_coeff_len: f32 = 40.0;
        let mapq_coeff_fac: f32 = mapq_coeff_len.ln();

        let mut mat = [0_i8; 25];
        fill_scmat(a, b, &mut mat);

        MemOpt {
            a,
            b,
            o_del: 6,
            e_del: 1,
            o_ins: 6,
            e_ins: 1,
            pen_unpaired: 17,
            pen_clip5: 5,
            pen_clip3: 5,
            w: 100,
            zdrop: 100,
            max_mem_intv: 20,
            t: 30,
            flag: 0,
            min_seed_len: 19,
            min_chain_weight: 0,
            max_chain_extend: 1 << 30,
            split_factor: 1.5,
            split_width: 10,
            max_occ: 500,
            max_chain_gap: 10_000,
            n_threads: 1,
            chunk_size: 10_000_000,
            mask_level: 0.50,
            drop_ratio: 0.50,
            xa_drop_ratio: 0.80,
            mask_level_redun: 0.95,
            mapq_coeff_len,
            mapq_coeff_fac,
            max_ins: 10_000,
            max_matesw: 50,
            max_xs_aux: 2,
            min_ins: 0,
            mat,
        }
    }
}

// ---------------------------------------------------------------------------
// Input sequence (bseq1_t in C)
// ---------------------------------------------------------------------------

/// A single input read, as produced by the FASTQ/FASTA reader.
#[derive(Debug, Clone, Default)]
pub struct BSeq {
    /// Read index within the current batch (0-based).
    pub id: i32,
    /// Read name (QNAME in SAM).
    pub name: String,
    /// Optional comment field from the FASTQ header.
    pub comment: String,
    /// ASCII nucleotide sequence (upper-case, e.g. `b"ACGT"`).
    pub seq: Vec<u8>,
    /// ASCII base-quality scores (Phred+33 encoded), or empty if FASTA input.
    pub qual: Vec<u8>,
    /// Accumulated SAM output line(s) for this read; may contain multiple records
    /// separated by `\n` for supplementary/secondary alignments.
    pub sam: String,
}

// ---------------------------------------------------------------------------
// Seed (mem_seed_t in C)
// ---------------------------------------------------------------------------

/// A single maximal exact match (MEM) seed, linking a query interval to a
/// reference position.
#[derive(Debug, Clone, Default)]
pub struct MemSeed {
    /// Reference begin position (0-based, in the concatenated forward strand).
    pub rbeg: i64,
    /// Query begin position (0-based).
    pub qbeg: i32,
    /// Seed length in bases.
    pub len: i32,
    /// Alignment score contributed by this seed.
    pub score: i32,
    /// Set to `true` once this seed has been processed during chaining/extension.
    pub done: bool,
}

// ---------------------------------------------------------------------------
// Chain (mem_chain_t in C)
// ---------------------------------------------------------------------------

/// A co-linear chain of seeds on the same reference strand, used as a
/// candidate alignment region before full Smith-Waterman extension.
#[derive(Debug, Clone, Default)]
pub struct MemChain {
    /// Index of the reference sequence this chain falls on.
    pub seqid: i32,
    /// Ordered list of seeds that form the chain.
    pub seeds: Vec<MemSeed>,
    /// Band width used (or estimated) for this chain.
    pub w: i32,
    /// Filtering status: 0 = filtered out, 1 = kept, 2 = already extended.
    pub kept: i32,
    /// Whether the chain maps to an ALT contig.
    pub is_alt: bool,
    /// Fraction of bases covered by repetitive seeds.
    pub frac_rep: f32,
    /// Total seed weight (coverage) for this chain; used during filtering.
    pub weight: i32,
}

// ---------------------------------------------------------------------------
// Alignment region (mem_alnreg_t in C)
// ---------------------------------------------------------------------------

/// A candidate alignment region produced from chain extension, before the
/// final CIGAR string is computed.
#[derive(Debug, Clone)]
pub struct MemAlnReg {
    /// Reference begin position of the aligned region (inclusive).
    pub rb: i64,
    /// Reference end position of the aligned region (exclusive).
    pub re: i64,
    /// Query begin position (inclusive).
    pub qb: i32,
    /// Query end position (exclusive).
    pub qe: i32,
    /// Reference sequence index (into `BntSeq::anns`).
    pub rid: i32,
    /// Best local SW score for this region.
    pub score: i32,
    /// True score corresponding to the aligned region (may be < `score` when
    /// the extension was clipped).
    pub truesc: i32,
    /// Second-best SW score (sub-optimal hit score).
    pub sub: i32,
    pub alt_sc: i32,
    /// SW score of the best tandem/chained hit.
    pub csub: i32,
    /// Approximate number of sub-optimal hits.
    pub sub_n: i32,
    /// Actual band width used during extension.
    pub w: i32,
    /// Number of reference bases covered by seeds in this region.
    pub seedcov: i32,
    /// Index of the primary hit that shadows this hit; negative if primary.
    pub secondary: i32,
    pub secondary_all: i32,
    /// Length of the seed that seeded this alignment.
    pub seedlen0: i32,
    /// Number of sub-alignments chained together to form this region.
    pub n_comp: i32,
    /// Whether this alignment is on an ALT contig.
    pub is_alt: bool,
    /// Fraction of the region covered by repetitive seeds.
    pub frac_rep: f32,
    /// Hash of the alignment region, used for deduplication.
    pub hash: u64,
}

impl Default for MemAlnReg {
    fn default() -> Self {
        MemAlnReg {
            rb: 0,
            re: 0,
            qb: 0,
            qe: 0,
            rid: 0,
            score: 0,
            truesc: 0,
            sub: 0,
            alt_sc: 0,
            csub: 0,
            sub_n: 0,
            w: 0,
            seedcov: 0,
            secondary: -1,
            secondary_all: -1,
            seedlen0: 0,
            n_comp: 0,
            is_alt: false,
            frac_rep: 0.0,
            hash: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// Finalized alignment (mem_aln_t in C)
// ---------------------------------------------------------------------------

/// A fully resolved alignment, including CIGAR string and SAM flag, as
/// returned by `mem_reg2aln`.
#[derive(Debug, Clone, Default)]
pub struct MemAln {
    /// Forward-strand mapping position (0-based).
    pub pos: i64,
    /// Reference sequence index; negative for unmapped reads.
    pub rid: i32,
    /// SAM FLAG field.
    pub flag: u16,
    /// `true` if the read is mapped to the reverse complement strand.
    pub is_rev: bool,
    /// `true` if the alignment is on an ALT contig.
    pub is_alt: bool,
    /// Mapping quality (MAPQ field in SAM).
    pub mapq: u8,
    /// Edit distance (NM tag).
    pub nm: i32,
    /// Number of CIGAR operations.
    pub n_cigar: i32,
    /// CIGAR encoded in BAM format: each element is `(op_length << 4) | op_code`
    /// where op codes follow the BAM convention (M=0, I=1, D=2, N=3, S=4, …).
    pub cigar: Vec<u32>,
    /// XA tag string listing alternative alignments, or empty.
    pub xa: String,
    /// Best SW score for this alignment.
    pub score: i32,
    /// Second-best SW score.
    pub sub: i32,
    pub alt_sc: i32,
}

// ---------------------------------------------------------------------------
// Paired-end insert-size statistics (mem_pestat_t in C)
// ---------------------------------------------------------------------------

/// Insert-size distribution statistics estimated from properly paired reads.
#[derive(Debug, Clone)]
pub struct MemPeStat {
    /// Lower bound of the "proper pair" insert-size range.
    pub low: i32,
    /// Upper bound of the "proper pair" insert-size range.
    pub high: i32,
    /// `true` when there were not enough concordant pairs to estimate the
    /// distribution reliably.
    pub failed: bool,
    /// Mean of the insert-size distribution.
    pub avg: f64,
    /// Standard deviation of the insert-size distribution.
    pub std: f64,
}

impl Default for MemPeStat {
    fn default() -> Self {
        MemPeStat {
            low: 0,
            high: 0,
            failed: true,
            avg: 0.0,
            std: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Super-Maximal Exact Match / FM-index interval (SMEM in C)
// ---------------------------------------------------------------------------

/// An FM-index BWT interval representing a super-maximal exact match (SMEM).
///
/// Corresponds to the `SMEM` struct (also called `bwtintv_t` in some bwa
/// versions).  The hit count is `l - k + 1`.
#[derive(Debug, Clone, Default)]
pub struct Smem {
    /// Reference sequence ID (not always populated; 0 when not set).
    pub rid: u32,
    /// Start position of the match in the query (inclusive).
    pub m: u32,
    /// End position of the match in the query (exclusive).
    pub n: u32,
    /// Lower bound of the BWT suffix-array interval (SA[k..=l] contain hits).
    pub k: i64,
    /// Upper bound of the BWT suffix-array interval.
    pub l: i64,
    /// Length of the exact match (number of bases), NOT the hit count.
    pub s: i64,
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mem_opt_default_fields() {
        let opt = MemOpt::default();
        assert_eq!(opt.a, 1);
        assert_eq!(opt.b, 4);
        assert_eq!(opt.o_del, 6);
        assert_eq!(opt.e_del, 1);
        assert_eq!(opt.o_ins, 6);
        assert_eq!(opt.e_ins, 1);
        assert_eq!(opt.pen_unpaired, 17);
        assert_eq!(opt.pen_clip5, 5);
        assert_eq!(opt.pen_clip3, 5);
        assert_eq!(opt.w, 100);
        assert_eq!(opt.zdrop, 100);
        assert_eq!(opt.t, 30);
        assert_eq!(opt.min_seed_len, 19);
        assert_eq!(opt.split_width, 10);
        assert_eq!(opt.max_occ, 500);
        assert_eq!(opt.max_chain_gap, 10_000);
        assert_eq!(opt.n_threads, 1);
        assert_eq!(opt.chunk_size, 10_000_000);
        assert!((opt.mask_level - 0.50).abs() < 1e-6);
        assert!((opt.drop_ratio - 0.50).abs() < 1e-6);
        assert!((opt.xa_drop_ratio - 0.80).abs() < 1e-6);
        assert!((opt.mask_level_redun - 0.95).abs() < 1e-6);
        assert_eq!(opt.max_ins, 10_000);
        assert_eq!(opt.max_matesw, 50);
        assert_eq!(opt.max_xs_aux, 2);
    }

    #[test]
    fn scoring_matrix_values() {
        let opt = MemOpt::default();
        // Diagonal entries should equal `a`.
        for i in 0..4usize {
            assert_eq!(opt.mat[i * 5 + i], opt.a as i8, "diagonal mismatch at {i}");
        }
        // Off-diagonal ACGT entries should equal `-b`.
        for i in 0..4usize {
            for j in 0..4usize {
                if i != j {
                    assert_eq!(
                        opt.mat[i * 5 + j],
                        -(opt.b as i8),
                        "off-diagonal mismatch at ({i},{j})"
                    );
                }
            }
        }
        // N row and N column should all be -1.
        for j in 0..5usize {
            assert_eq!(opt.mat[4 * 5 + j], -1, "N row mismatch at col {j}");
            assert_eq!(opt.mat[j * 5 + 4], -1, "N col mismatch at row {j}");
        }
    }

    #[test]
    fn aln_reg_default_secondary() {
        let reg = MemAlnReg::default();
        assert_eq!(reg.secondary, -1);
        assert_eq!(reg.secondary_all, -1);
    }

    #[test]
    fn pe_stat_default_failed() {
        let stat = MemPeStat::default();
        assert!(stat.failed);
        assert_eq!(stat.low, 0);
        assert_eq!(stat.high, 0);
    }
}
