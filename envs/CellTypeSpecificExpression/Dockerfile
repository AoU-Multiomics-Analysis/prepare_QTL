FROM mambaorg/micromamba:2.9.0-ubuntu22.04

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/CellTypeSpecificExpression/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/cell_type_specific_expression /opt/prepare_qtl/scripts/cell_type_specific_expression
ENV PATH=/opt/conda/bin:$PATH
ENV R_LIBS_USER=/opt/conda/lib/R/library
ENV CELL_TYPE_SPECIFIC_EXPRESSION_ROOT=/opt/prepare_qtl/scripts/cell_type_specific_expression

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript -e 'stopifnot(as.character(packageVersion("dtangle")) == "2.0.10")' && \
    Rscript -e 'stopifnot(as.character(packageVersion("TCA")) == "1.2.1")'

WORKDIR /workspace
