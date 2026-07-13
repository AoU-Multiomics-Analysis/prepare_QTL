# Slim runtime for the high-throughput, Rust-only methylation tasks. It avoids
# the R/Bioconductor environment used by the downstream cohort analyses.
FROM rust:1.85-slim-bullseye AS builder

WORKDIR /build/rust/methylation_filter
COPY rust/methylation_filter/Cargo.toml rust/methylation_filter/Cargo.lock ./
COPY rust/methylation_filter/src ./src
RUN cargo test --release && cargo build --release

WORKDIR /build/rust/methylation_merge
COPY rust/methylation_merge/Cargo.toml rust/methylation_merge/Cargo.lock ./
COPY rust/methylation_merge/src ./src
RUN cargo test --release && cargo build --release

FROM debian:bullseye-slim
COPY --from=builder /build/rust/methylation_filter/target/release/methylation-filter /usr/local/bin/methylation-filter
COPY --from=builder /build/rust/methylation_merge/target/release/methylation-chromosome-merge /usr/local/bin/methylation-chromosome-merge
CMD ["bash"]
