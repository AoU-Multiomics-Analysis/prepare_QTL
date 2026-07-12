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
COPY scripts/ /opt/prepare_qtl/scripts/
## Default command to run when the container starts
CMD ["bash"]
