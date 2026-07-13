# prepare_QTL

Pipeline for preparing molecular phenotype and genotype data for QTL (Quantitative Trait Loci) analysis. Supports expression (eQTL), splicing (sQTL), and proteomics (pQTL) data types.

## Repository Structure

```
prepare_QTL/
├── scripts/
│   ├── common/     # Shared phenotype-PC, covariate, and residualization utilities
│   ├── expression/ # eQTL phenotype preparation
│   ├── methylation/ # PacBio 5mC filtering, aggregation, annotation, and correlation QC
│   ├── proteomics/ # Olink normalization and pQTL preparation
│   └── splicing/   # sQTL phenotype preparation
├── workflows/      # WDL workflows for running analyses on a cloud platform
├── envs/           # Docker environments
└── docs/           # Detailed documentation
```

## Documentation

- [R scripts](docs/scripts.md): Command-line script inputs, outputs, and processing behavior.
- [Phenotype normalization and filtering](docs/phenotype-normalization-filtering.md): Modality-specific filtering, normalization, transformations, and WGCNA outlier removal.
- [Molecular QTL workflows](docs/molecular-qtl-workflows.md): eQTL, sQTL, pQTL, methylation, proteomics normalization, phenotype PC, and covariate merge WDLs.
- [PacBio 5mC QTL workflow](docs/methylation-qtl.md): pb-CpG-tools inputs, QC, sharding, site metadata, and TensorQTL phenotype output.
- [Genotype workflows](docs/genotype-workflows.md): VCF, PLINK, genotype PC, allele frequency, and dosage WDLs.
- [Docker environment](docs/docker.md): Docker image location and included R package dependencies.

## Main Workflows

- [`workflows/expression/prepare_eQTL.wdl`](workflows/expression/prepare_eQTL.wdl): Prepares INT, scaled, and raw expression BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/splicing/prepare_sQTL.wdl`](workflows/splicing/prepare_sQTL.wdl): Prepares INT, scaled, and raw splice BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/proteomics/prepare_pQTL.wdl`](workflows/proteomics/prepare_pQTL.wdl): Prepares INT, scaled, and raw proteomics BED files, plus phenotype PCs, optional QTL covariates, and optional residualized BEDs for INT and scaled outputs.
- [`workflows/proteomics/normalize_pQTL.wdl`](workflows/proteomics/normalize_pQTL.wdl): Median-normalizes Olink NPX parquet files before pQTL preparation.
- [`workflows/methylation/ProcessMethylationSample.wdl`](workflows/methylation/ProcessMethylationSample.wdl) and [`workflows/methylation/AggregateMethylationCohort.wdl`](workflows/methylation/AggregateMethylationCohort.wdl): Terra-table sample and cohort entry points for pb-CpG-tools 5mC QTL preparation. Cohort aggregation accepts one compact output manifest rather than 22 Terra file arrays. [`workflows/methylation/merge_methylation.wdl`](workflows/methylation/merge_methylation.wdl) remains the source-BED manifest/shard wrapper.
- [`workflows/genotype/prepare_VCF.wdl`](workflows/genotype/prepare_VCF.wdl): Prepares genotype data from an All of Us Hail MatrixTable.

See the [workflow catalog](workflows/README.md) for public entry points, internal building blocks, and maintenance conventions.

## Common Options

The prepare scripts and workflows for eQTL, pQTL, and sQTL share this output pattern:

- `.INT`: Rank-based inverse normal transformed molecular phenotypes.
- `.scaled`: Centered and scaled molecular phenotypes. Expression CPMs are transformed with `log2(CPM + 1)` before centering/scaling; proteomics and splicing values are centered/scaled directly.
- `.raw`: Untransformed phenotype values after sample/feature filtering and BED formatting. Raw BEDs are emitted as workflow outputs but are not used for phenotype PCs or covariate merging.
- Connectivity outliers: `.INT` and `.scaled` BEDs have WGCNA sample connectivity outliers removed before downstream PC, covariate, or residualization steps. Removed samples are written to `*.connectivity_outliers.tsv`; raw BEDs keep all samples after the initial sample-list filter.
- `AdditionalCovariates`: Optional WDL input for eQTL, pQTL, and sQTL prepare workflows. When provided, the workflow runs [`workflows/common/MergeCovariates.wdl`](workflows/common/MergeCovariates.wdl) for both `.INT` and `.scaled` phenotype PCs.
- `ResidualizeNormalizedInputs`: Optional WDL toggle for eQTL, pQTL, and sQTL prepare workflows. When `true`, the workflow writes residualized BEDs for `.INT` and `.scaled` inputs; raw BEDs are not residualized.

See the linked docs above for full input and output details.
