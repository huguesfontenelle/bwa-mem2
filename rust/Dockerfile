# ─────────────────────────────────────────────
# Stage 1: build bwa-mem2 (Rust)
# ─────────────────────────────────────────────
FROM rust:1.85-slim AS builder

WORKDIR /build

# Cache dependency compilation separately from source
COPY Cargo.toml Cargo.lock ./
# Create a dummy main so cargo can compile deps
RUN mkdir src && echo 'fn main(){}' > src/main.rs
RUN cargo build --release
RUN rm src/main.rs

# Now copy real source and do the final build
COPY src/ src/
# Touch main.rs so cargo knows it changed
RUN touch src/main.rs
RUN cargo build --release

# ─────────────────────────────────────────────
# Stage 2: build samtools 1.20 from source
# ─────────────────────────────────────────────
FROM debian:bookworm-slim AS samtools-builder

ARG SAMTOOLS_VERSION=1.20

RUN apt-get update && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        bzip2 \
        ca-certificates \
        curl \
        gcc \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        make \
        perl \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build/samtools

RUN curl -fsSL \
    "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    | tar -xj --strip-components=1

RUN ./configure --prefix=/opt/samtools --without-curses \
    && make -j"$(nproc)" \
    && make install

# ─────────────────────────────────────────────
# Stage 3: minimal runtime image
# ─────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        libbz2-1.0 \
        libcurl4 \
        liblzma5 \
        libncurses6 \
        libssl3 \
        procps \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder      /build/target/release/bwa-mem2 /usr/local/bin/bwa-mem2
COPY --from=samtools-builder /opt/samtools/bin/samtools  /usr/local/bin/samtools
