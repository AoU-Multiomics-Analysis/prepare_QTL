FROM mambaorg/micromamba:2.9.0-ubuntu22.04

LABEL org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/prepare_QTL"
LABEL org.opencontainers.image.description="Batched RNA-SeQC table aggregation for Terra"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/RNASeQCAggregation/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/expression/merge_rnaseqc.py /opt/prepare_qtl/scripts/expression/merge_rnaseqc.py
# Standalone gsutil needs this setting to request credentials from VM metadata.
COPY envs/RNASeQCAggregation/boto.cfg /etc/boto.cfg
# Terra can bypass the micromamba entrypoint. These tools need no activation hook.
ENV PATH=/opt/conda/bin:$PATH

WORKDIR /tmp
CMD ["bash"]
