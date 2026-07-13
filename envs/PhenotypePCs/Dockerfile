## Build the streaming per-sample methylation filter as a compiled workflow tool.
## Bullseye keeps the resulting binary compatible with the micromamba runtime.
FROM rust:1.85-slim-bullseye AS methylation_filter_builder
WORKDIR /build/rust/methylation_filter
COPY rust/methylation_filter/Cargo.toml ./
COPY rust/methylation_filter/src ./src
RUN cargo test --release && cargo build --release

## Use micromamba as the base image
FROM mambaorg/micromamba:1.5.3

## Set up environment
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
## Create a new environment and install packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::r-data.table \
    conda-forge::r-arrow \
    conda-forge::r-olinkanalyze \
    conda-forge::r-optparse \
    conda-forge::google-cloud-sdk \
    conda-forge::crcmod \
    conda-forge::r-janitor \
    conda-forge::r-ggupset \
    bioconda::bioconductor-pcatools \
    conda-forge::r-rnomni \
    bioconda::r-wgcna \
    bioconda::bioconductor-edger \
    bioconda::bioconductor-biomart \
    bioconda::r-biomartr \
    bioconda::bioconductor-rtracklayer \
    bioconda::bioconductor-plyranges \
    conda-forge::r-patchwork \
    bioconda::htslib
RUN gsutil version -l | grep -q 'compiled crcmod: True'
COPY --from=methylation_filter_builder /build/rust/methylation_filter/target/release/methylation-filter /usr/local/bin/methylation-filter
COPY scripts/ /opt/prepare_qtl/scripts/
# Keep the legacy flat location while Terra workflows transition between the
# mutable :main image tag and a newly published image. Source files remain
# organized under /opt/prepare_qtl/scripts/.
COPY scripts/common/*.R scripts/expression/*.R scripts/methylation/*.R scripts/proteomics/*.R scripts/splicing/*.R /tmp/
## Default command to run when the container starts
CMD ["bash"]
