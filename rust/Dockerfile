# ─────────────────────────────────────────────
# Stage 1: build
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
# Stage 2: minimal runtime image
# ─────────────────────────────────────────────
FROM debian:bookworm-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/bwa-mem2 /usr/local/bin/bwa-mem2

ENTRYPOINT ["bwa-mem2"]
CMD ["--help"]
