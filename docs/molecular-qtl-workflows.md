# Molecular QTL Workflows

[Back to main README](../README.md)

These WDL workflows prepare molecular phenotype data for expression, splicing, and proteomics QTL analyses. They are designed to run on a cloud platform such as Terra and wrap the R scripts documented in [R scripts](scripts.md).

## Shared Prepare Workflow Outputs

The eQTL, pQTL, and sQTL prepare workflows now compute both molecular phenotype transformations:

- `.INT`: Rank-based inverse normal transformed molecular phenotypes.
- `.scaled`: Centered and scaled molecular phenotypes.
- `.raw`: Untransformed phenotype values after sample/feature filtering and BED formatting.

Each workflow computes phenotype PCs separately for the `.INT` and `.scaled` outputs only. Raw BED files are emitted as workflow outputs but are not used for phenotype PCs or covariate merging. `AdditionalCovariates` is an optional TSV of covariates with a `sample_id` column. When provided, the workflow runs [`MergeCovariates.wdl`](../workflows/MergeCovariates.wdl) twice to merge those covariates with the `.INT` and `.scaled` phenotype PCs.

## `workflows/prepare_eQTL.wdl`

End-to-end workflow for preparing gene expression data for eQTL analysis.

**Steps:**
1. Runs `PrepareExpression.R` to produce `.INT`, `.scaled`, and `.raw` expression BED files.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` expression BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.

**Inputs:** Raw count GCT file, GENCODE GTF, sample list, output prefix, optional additional covariates TSV, resource parameters.

**Outputs:** `.expression.INT.bed.gz`, `.expression.scaled.bed.gz`, `.expression.raw.bed.gz`, phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, and optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv`.

## `workflows/prepare_pQTL.wdl`

End-to-end workflow for preparing Olink proteomics data for pQTL analysis.

**Steps:**
1. Runs `PrepareProteomics.R` to produce `.INT`, `.scaled`, and `.raw` protein BED files.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` protein BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.

**Inputs:** Olink proteomics data file, GENCODE GTF, sample list, output prefix, optional additional covariates TSV, resource parameters.

**Outputs:** `.protein.INT.bed.gz`, `.protein.scaled.bed.gz`, `.protein.raw.bed.gz`, phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, and optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv`.

## `workflows/normalize_pQTL.wdl`

Workflow that median-normalizes Olink NPX parquet files before pQTL preparation. This workflow is registered in `.dockstore.yml` as `normalize_pQTL`.

**Inputs:** `Array[File] OlinkData` containing Olink NPX parquet files, output prefix, reference plate ID, resource parameters.

**Outputs:** Median-normalized Olink TSV and filtered long-format proteomics TSV. The filtered output can be passed directly to [`workflows/prepare_pQTL.wdl`](../workflows/prepare_pQTL.wdl) as `ProteomicData`.

## `workflows/prepare_sQTL.wdl`

End-to-end workflow for preparing splice junction data for sQTL analysis.

**Steps:**
1. Runs `PrepareSpliceData.R` to produce `.INT`, `.scaled`, and `.raw` splice BED files.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` splice BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.

**Inputs:** LeafCutter BED file, sample list, output prefix, optional additional covariates TSV, resource parameters.

**Outputs:** `.splicing.INT.bed.gz`, `.splicing.scaled.bed.gz`, `.splicing.raw.bed.gz`, phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, and optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv`.

## `workflows/calculate_phenotypePCs.wdl`

Workflow that computes phenotype PCs from any normalized molecular phenotype BED file.

**Steps:**
1. Runs `calculate_PCs.R` on the input BED file using the Gavish-Donoho method to select the number of PCs.

**Inputs:** Normalized BED file, output prefix, optional output suffix, resource parameters.

**Outputs:** Phenotype PCs TSV (`<OutputPrefix>_phenotype_PCs<OutputSuffix>.tsv`).

## `workflows/MergeCovariates.wdl`

Workflow that merges additional covariates, such as genotype PCs, and molecular phenotype PCs into a single covariate file ready for tensorQTL.

**Steps:**
1. Runs `MergeCovariates.R` to inner-join additional covariates and molecular PCs, then transpose the result.

**Inputs:** Additional covariates TSV with `sample_id`, molecular PCs TSV, output prefix, optional output suffix.

**Outputs:** Combined QTL covariate file (`<OutputPrefix>_QTL_covariates<OutputSuffix>.tsv`).
