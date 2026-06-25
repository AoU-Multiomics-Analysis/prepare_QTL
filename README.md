# prepare_QTL

Pipeline for preparing molecular phenotype and genotype data for QTL (Quantitative Trait Loci) analysis. Supports expression (eQTL), splicing (sQTL), and proteomics (pQTL) data types.

## Repository Structure

```
prepare_QTL/
├── scripts/        # R scripts for data processing and normalization
├── workflows/      # WDL workflows for running analyses on a cloud platform
├── envs/           # Docker environments
└── docs/           # Detailed documentation
```

## Documentation

- [R scripts](docs/scripts.md): Command-line script inputs, outputs, and processing behavior.
- [Molecular QTL workflows](docs/molecular-qtl-workflows.md): eQTL, sQTL, pQTL, proteomics normalization, phenotype PC, and covariate merge WDLs.
- [Genotype workflows](docs/genotype-workflows.md): VCF, PLINK, genotype PC, allele frequency, and dosage WDLs.
- [Docker environment](docs/docker.md): Docker image location and included R package dependencies.

## Main Workflows

- [`workflows/prepare_eQTL.wdl`](workflows/prepare_eQTL.wdl): Prepares INT, scaled, and raw expression BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/prepare_sQTL.wdl`](workflows/prepare_sQTL.wdl): Prepares INT, scaled, and raw splice BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/prepare_pQTL.wdl`](workflows/prepare_pQTL.wdl): Prepares INT, scaled, and raw proteomics BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/normalize_pQTL.wdl`](workflows/normalize_pQTL.wdl): Median-normalizes Olink NPX parquet files before pQTL preparation.
- [`workflows/prepare_VCF.wdl`](workflows/prepare_VCF.wdl): Prepares genotype data from an All of Us Hail MatrixTable.

## Common Options

The prepare scripts and workflows for eQTL, pQTL, and sQTL share this output pattern:

- `.INT`: Rank-based inverse normal transformed molecular phenotypes.
- `.scaled`: Centered and scaled molecular phenotypes.
- `.raw`: Untransformed phenotype values after sample/feature filtering and BED formatting. Raw BEDs are emitted as workflow outputs but are not used for phenotype PCs or covariate merging.
- `AdditionalCovariates`: Optional WDL input for eQTL, pQTL, and sQTL prepare workflows. When provided, the workflow runs `MergeCovariates.wdl` for both `.INT` and `.scaled` phenotype PCs.
- `ResidualizeNormalizedInputs`: Optional WDL toggle for eQTL, pQTL, and sQTL prepare workflows. When `true`, the workflow writes residualized BEDs for `.INT` and `.scaled` inputs; raw BEDs are not residualized.

See the linked docs above for full input and output details.
