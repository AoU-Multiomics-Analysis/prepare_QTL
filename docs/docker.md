# Docker Environment

[Back to main README](../README.md)

The Docker image used by most WDL tasks is defined in [`envs/PhenotypePCs/Dockerfile`](../envs/PhenotypePCs/Dockerfile) and published as:

```text
ghcr.io/aou-multiomics-analysis/prepare_qtl:main
```

The image is built automatically on every push or pull request to `main` through the GitHub Actions workflow in [`.github/workflows/docker-image.yml`](../.github/workflows/docker-image.yml).

## Included R Packages

The image includes the following R packages used by the scripts and WDL tasks:

- `tidyverse`, `data.table`, `arrow`, `OlinkAnalyze`, `optparse`, `janitor`, `WGCNA` (`bioconda::r-wgcna`)
- `PCAtools`, `RNOmni`, `edgeR`
- `biomaRt`, `biomaRtr`, `rtracklayer`, `plyranges`
- `patchwork`
