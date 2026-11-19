# syntax=docker/dockerfile:1.7
FROM mambaorg/micromamba:1.5.3

ENV PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin
RUN micromamba install --yes --name base --strict-channel-priority \
      --channel conda-forge --channel bioconda \
      python=3.11.13 htslib=1.22.1 \
    && micromamba clean --all --yes

# PacBio's Linux x86_64 release, verified before extraction.
ADD --checksum=sha256:442e91524236d42d139d4e129ed643076d5fafa545d01b65d34a2c53cf3c3683 --chmod=0644 \
    https://github.com/PacificBiosciences/MethBat/releases/download/v1.2.0/methbat-v1.2.0-x86_64-unknown-linux-gnu.tar.gz \
    /tmp/methbat.tar.gz
RUN tar -xzf /tmp/methbat.tar.gz --strip-components=1 -C /opt/conda/bin \
      methbat-v1.2.0-x86_64-unknown-linux-gnu/methbat \
    && test "$(methbat --version)" = "methbat 1.2.0-7a582da"

COPY --chmod=0755 methbat_prepare_pileups.py methbat_compare_methods.py /opt/conda/bin/
CMD ["bash"]
