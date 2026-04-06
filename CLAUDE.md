# CLAUDE.md

## Project Overview

bwa-mem2 is a high-performance DNA sequence aligner, the next version of bwa-mem. It produces alignment identical to bwa and is ~1.3-3.1x faster. Written in C/C++ with SIMD optimizations (SSE2, AVX2, AVX-512).

## Branch: `rust`

This branch is for porting or wrapping bwa-mem2 functionality in Rust. The Rust port lives in `rust/`.

---

## Rust Port

### Requirements

- Rust 1.70+ (`rustup` recommended)
- Cargo (included with Rust)

### Build

```sh
cd rust

# Debug build
cargo build

# Release build (optimised)
cargo build --release

# Run tests
cargo test
```

The compiled binary is at:

```
rust/target/release/bwa-mem2       # release
rust/target/debug/bwa-mem2         # debug
```

### Usage

#### 1. Build an index

```sh
bwa-mem2 index ref.fa
# creates ref.fa.pac, ref.fa.ann, ref.fa.amb, ref.fa.bwt.2bit.64

# Custom index prefix
bwa-mem2 index -p myref ref.fa
```

#### 2. Align reads (single-end)

```sh
bwa-mem2 mem ref.fa reads.fq > out.sam
bwa-mem2 mem ref.fa reads.fq.gz > out.sam   # gzip input
```

#### 3. Align reads (paired-end)

```sh
bwa-mem2 mem ref.fa read1.fq read2.fq > out.sam
```

#### 4. Common options

| Flag | Default | Description |
|------|---------|-------------|
| `-t INT` | 1 | Number of threads |
| `-k INT` | 19 | Minimum seed length |
| `-w INT` | 100 | Band width |
| `-d INT` | 100 | Off-diagonal Z-drop |
| `-c INT` | 500 | Skip seeds with >INT occurrences |
| `-A INT` | 1 | Match score |
| `-B INT` | 4 | Mismatch penalty |
| `-O INT` | 6 | Gap open penalty |
| `-E INT` | 1 | Gap extension penalty |
| `-L INT` | 5 | Clipping penalty |
| `-U INT` | 17 | Unpaired penalty (PE mode) |
| `-R STR` | — | Read group line (e.g. `@RG\tID:foo\tSM:bar`) |
| `-o FILE` | stdout | Output SAM file |
| `-a` | off | Output all alignments |
| `-M` | off | Mark shorter hits as secondary |
| `-Y` | off | Use soft clipping for supplementary |

#### 5. Print version

```sh
bwa-mem2 version
```

### Rust Source Layout

```
rust/
├── Cargo.toml
└── src/
    ├── main.rs          — CLI entry point (clap), batch pipeline, SAM header
    ├── types.rs         — Core structs: MemOpt, BSeq, MemSeed, MemChain,
    │                      MemAlnReg, MemAln, MemPeStat, Smem
    ├── bntseq.rs        — Reference genome packing, FASTA→PAC, .ann/.amb I/O
    ├── bwt.rs           — FM-index (CpOcc, FmIndex), occ/LF-mapping, I/O
    ├── kseq.rs          — Streaming FASTA/FASTQ reader (plain + gzip)
    ├── index.rs         — SA-IS suffix array + FM-index construction
    ├── fmi_search.rs    — Bi-directional SMEM finding via FM-index
    ├── sw.rs            — Banded Smith-Waterman (local + global + CIGAR)
    ├── mem.rs           — Seeding, chaining, filtering, extension, MAPQ
    ├── sam.rs           — SAM output, NM/MD tags, cigar_to_str, rev_comp
    └── pe.rs            — Paired-end insert-size estimation, pair scoring
```

---

## C/C++ Build

```sh
# Default build
make

# Multi-architecture build (native, SSE4.1, AVX2, AVX-512)
make multi

# With Intel compiler
make CXX=icpc multi
```

## Repository Structure

```
rust/         Rust port (this branch)
src/          C/C++ source files
ext/          External dependencies (safestringlib, etc.)
test/         Test data and scripts
images/       Documentation images
```

## Key C/C++ Source Files

- `src/fastmap.cpp` — main entry point, index and mem subcommands
- `src/bwamem.cpp` — core mem algorithm
- `src/FMI_search.cpp` — FM-index search
- `src/bandedSWA.cpp` — banded Smith-Waterman alignment
- `src/bwtindex.cpp` — index construction

## Git Submodules

```sh
git submodule init
git submodule update
```

## Development Notes

- All changes on this branch should be committed and pushed to `origin/rust`
- Rust code goes in `rust/`; use `cargo` for build
- Parallelism via `rayon`; no unsafe code except the PE reader lifetime workaround
