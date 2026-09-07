FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --strict-channel-priority --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/ /opt/trans-ld-regions/scripts/
ENV PATH=/opt/conda/bin:$PATH
ENV PYTHONPATH=/opt/trans-ld-regions/scripts
ENV PYTHONUNBUFFERED=1
CMD ["python3", "-m", "trans_ld_regions.prepare", "--help"]
