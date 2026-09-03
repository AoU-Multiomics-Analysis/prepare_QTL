FROM mambaorg/micromamba:2.9.0-ubuntu22.04

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/CellTypeSpecificExpression/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/cell_type_specific_expression /opt/prepare_qtl/scripts/cell_type_specific_expression
ENV PATH=/opt/conda/bin:$PATH
ENV R_LIBS_USER=/opt/conda/lib/R/library
ENV CELL_TYPE_SPECIFIC_EXPRESSION_ROOT=/opt/prepare_qtl/scripts/cell_type_specific_expression

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# HSPE is not on CRAN. Verify the pinned upstream source before installation.
RUN Rscript -e 'archive <- tempfile(fileext = ".tar.gz"); \
    download.file("https://gjhunt.github.io/hspe/hspe_0.1.tar.gz", archive, mode = "wb"); \
    stopifnot(digest::digest(file = archive, algo = "sha256", serialize = FALSE) == \
      "f6c7c27ba2aec9b042cb8f793d7b7565023461812cdebdc0228b240d67b68c78"); \
    install.packages(archive, repos = NULL, type = "source"); \
    unlink(archive); \
    stopifnot(as.character(packageVersion("hspe")) == "0.1"); \
    stopifnot(as.character(packageVersion("DEoptimR")) == "1.1.4"); \
    stopifnot(as.character(packageVersion("nloptr")) == "2.2.1")' && \
    Rscript -e 'stopifnot(as.character(packageVersion("TCA")) == "1.2.1")'

WORKDIR /workspace
