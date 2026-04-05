# CLAUDE.md

## Project Overview

bwa-mem2 is a high-performance DNA sequence aligner, the next version of bwa-mem. It produces alignment identical to bwa and is ~1.3-3.1x faster. Written in C/C++ with SIMD optimizations (SSE2, AVX2, AVX-512).

## Branch: `rust`

This branch is for porting or wrapping bwa-mem2 functionality in Rust.

## Build

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
src/          C/C++ source files
ext/          External dependencies (safestringlib, etc.)
test/         Test data and scripts
images/       Documentation images
```

## Key Source Files

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
- Rust code goes in a new `rust/` directory at the repo root
- Use `cargo` for Rust build; interop with C via FFI as needed
