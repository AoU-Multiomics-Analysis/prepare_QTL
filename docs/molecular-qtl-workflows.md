# Molecular QTL Workflows

[Back to main README](../README.md)

These WDL workflows prepare molecular phenotype data for expression, splicing, and proteomics QTL analyses. They are designed to run on a cloud platform such as Terra and wrap the R scripts documented in [R scripts](scripts.md). For modality-specific normalization and filtering details, see [Phenotype normalization and filtering](phenotype-normalization-filtering.md).

## Shared Prepare Workflow Outputs

For individual RNA-SeQC outputs, first use [`rnaseqc2_aggregate_batched.wdl`](../workflows/expression/rnaseqc2_aggregate_batched.wdl) to build cohort GCT and metrics files from one sample manifest. See the [RNA-SeQC aggregation guide](rnaseqc2-aggregation.md) for file naming, batch settings, and validation.

The eQTL, pQTL, and sQTL prepare workflows now compute both molecular phenotype transformations:

- `.INT`: Rank-based inverse normal transformed molecular phenotypes.
- `.scaled`: Centered and scaled molecular phenotypes. Raw-count expression mode transforms CPMs with `log2(CPM + 1)` before centering and scaling. Pre-normalized log2-CPM expression, proteomics, and splicing values are centered and scaled directly.
- `.raw`: Untransformed phenotype values after sample/feature filtering and BED formatting.

Each workflow computes phenotype PCs separately for the `.INT` and `.scaled` outputs only. For each transformed branch, the workflow emits both the existing Gavish-Donoho-selected phenotype-PC TSV and a matching `.all.tsv` file containing every available rotated PC from the same PCA run. Raw BED files are emitted as workflow outputs but are not used for phenotype PCs, covariate merging, or residualization. `AdditionalCovariates` is an optional TSV of covariates with a `sample_id` column. When provided, the workflow runs [`MergeCovariates.wdl`](../workflows/common/MergeCovariates.wdl) twice to merge those covariates with the selected `.INT` and `.scaled` phenotype-PC TSVs.

The prepare scripts remove WGCNA sample connectivity outliers from `.INT` and `.scaled` matrices before those BED files are emitted. For each transformed branch, the scripts compute a sample-sample biweight midcorrelation matrix with `WGCNA::bicor(..., use = "pairwise.complete.obs")`, transform it to adjacency with `0.5 + 0.5 * correlation`, calculate sample connectivity with `WGCNA::fundamentalNetworkConcepts()`, and remove samples with connectivity `Z_score < -3`. Removed samples are written as `connectivity_outliers.tsv` outputs with `SampleID` and `Z_score` columns. If there are fewer than 3 samples, fewer than 2 features, or zero/undefined connectivity variance, the scripts keep all samples and write an empty outlier TSV. Raw BED files keep all samples after the initial sample-list filter.

Set `ResidualizeNormalizedInputs` to `true` to run [`ResidualizePhenotypes.wdl`](../workflows/common/ResidualizePhenotypes.wdl) for the `.INT` and `.scaled` BED files. When merged covariates are available, the residualization task regresses each phenotype row on the corresponding merged covariates and then centers/scales the residuals. Without merged covariates, the task only centers/scales the input phenotype rows.

## `workflows/expression/prepare_eQTL.wdl`

End-to-end workflow for preparing gene expression data for eQTL analysis.

**Steps:**
1. Runs `PrepareExpression.R` to produce `.INT`, `.scaled`, and `.raw` expression BED files. Raw-count mode applies TMM and CPM; its `.scaled` branch applies `log2(CPM + 1)` before centering/scaling. Alternatively, `Log2CpmBed` accepts pre-normalized log2 CPM values, skips count and CPM processing, rank-normalizes those values for `.INT`, and centers/scales them directly for `.scaled`.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/common/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` expression BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/common/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.
4. Optionally runs [`ResidualizePhenotypes.wdl`](../workflows/common/ResidualizePhenotypes.wdl) for the `.INT` and `.scaled` BED files when `ResidualizeNormalizedInputs` is `true`.

**Inputs:** Either a raw-count GCT plus GENCODE GTF or a pre-normalized log2-CPM BED, sample list, output prefix, optional additional covariates TSV, residualization toggle, resource parameters.

**Outputs:** `.expression.INT.bed.gz`, `.expression.scaled.bed.gz`, `.expression.raw.bed.gz`, connectivity outlier TSVs for `.INT` and `.scaled`, selected phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, full phenotype-PC matrices ending in `.INT.all.tsv` and `.scaled.all.tsv`, optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv` that continue to use the selected PC files, and optionally residualized BEDs ending in `.residualized.bed.gz`.

## `workflows/proteomics/prepare_pQTL.wdl`

End-to-end workflow for preparing Olink proteomics data for pQTL analysis.

**Steps:**
1. Runs `PrepareProteomics.R` to produce `.INT`, `.scaled`, and `.raw` protein BED files.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/common/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` protein BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/common/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.
4. Optionally runs [`ResidualizePhenotypes.wdl`](../workflows/common/ResidualizePhenotypes.wdl) for the `.INT` and `.scaled` BED files when `ResidualizeNormalizedInputs` is `true`.

**Inputs:** Olink proteomics data file, GENCODE GTF, sample list, output prefix, optional additional covariates TSV, residualization toggle, resource parameters.

**Outputs:** `.protein.INT.bed.gz`, `.protein.scaled.bed.gz`, `.protein.raw.bed.gz`, connectivity outlier TSVs for `.INT` and `.scaled`, selected phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, full phenotype-PC matrices ending in `.INT.all.tsv` and `.scaled.all.tsv`, optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv` that continue to use the selected PC files, and optionally residualized BEDs ending in `.residualized.bed.gz`.

## `workflows/proteomics/normalize_pQTL.wdl`

Workflow that median-normalizes Olink NPX parquet files before pQTL preparation. This workflow is registered in `.dockstore.yml` as `normalize_pQTL`.

**Inputs:** `Array[File] OlinkData` containing Olink NPX parquet files, output prefix, reference plate ID, resource parameters.

**Outputs:** Median-normalized Olink TSV and filtered long-format proteomics TSV. The filtered output can be passed directly to [`workflows/proteomics/prepare_pQTL.wdl`](../workflows/proteomics/prepare_pQTL.wdl) as `ProteomicData`.

## Methylation workflows

`workflows/methylation/prepare_mQTL.wdl` is the entry point for an existing merged methylation BED. It accepts a headerless one-column sample list or a list headed by `sample_id`, `SampleID`, or `ID`; every requested sample must exist in the BED. It removes `.meth_region_stats` suffixes, keeps autosomal features below the missingness threshold, and mean-imputes retained missing values. At least two samples and two retained features are required. It writes raw, INT, and scaled BED files. Connectivity filtering applies independently to the INT and scaled files, while the raw file keeps all selected samples. The workflow calculates PCs for the two normalized files and optionally merges them with additional covariates.

`workflows/methylation/ProcessMethylationSample.wdl` runs per sample from direct Terra table fields (`SampleID` and `MethylationBed`), with no external manifest. It writes sample QC plus 22 autosome-specific QC-flagged call tables.

`workflows/methylation/AggregateMethylationCohort.wdl` consumes one compact cohort manifest containing those per-sample outputs, reconstructs the cohort sample list from the QC files, evaluates sample-presence and methylation-MAD filters in parallel per autosome, and aggregates final metadata and BEDs. Missing values among retained sites are imputed with the cohort feature mean before it writes raw beta-value and inverse-normalized TensorQTL BEDs, calculates both selected and full phenotype-PC outputs from the INT BED, and optionally merges the selected phenotype PCs with additional covariates.

`workflows/methylation/merge_methylation.wdl` is the source-BED manifest/shard wrapper for batch processing. It retains concurrent `gsutil` localization and uses an internal array-based cohort implementation, avoiding Terra's external workflow-input size limit.

See the [PacBio 5mC QTL workflow guide](methylation-qtl.md) for table bindings, manifest input, configuration, parallel filtering behavior, QC logic, and outputs.

## `workflows/splicing/prepare_sQTL.wdl`

End-to-end workflow for preparing splice junction data for sQTL analysis.

**Steps:**
1. Runs `PrepareSpliceData.R` to produce `.INT`, `.scaled`, and `.raw` splice BED files.
2. Runs `calculate_PCs.R` through [`calculate_phenotypePCs.wdl`](../workflows/common/calculate_phenotypePCs.wdl) separately for the `.INT` and `.scaled` splice BED files.
3. Optionally runs [`MergeCovariates.wdl`](../workflows/common/MergeCovariates.wdl) separately for the `.INT` and `.scaled` phenotype PCs when `AdditionalCovariates` is provided.
4. Optionally runs [`ResidualizePhenotypes.wdl`](../workflows/common/ResidualizePhenotypes.wdl) for the `.INT` and `.scaled` BED files when `ResidualizeNormalizedInputs` is `true`.

**Inputs:** LeafCutter BED file, sample list, output prefix, optional additional covariates TSV, residualization toggle, resource parameters.

**Outputs:** `.splicing.INT.bed.gz`, `.splicing.scaled.bed.gz`, `.splicing.raw.bed.gz`, connectivity outlier TSVs for `.INT` and `.scaled`, selected phenotype PCs ending in `.INT.tsv` and `.scaled.tsv`, full phenotype-PC matrices ending in `.INT.all.tsv` and `.scaled.all.tsv`, optionally merged QTL covariates ending in `.INT.tsv` and `.scaled.tsv` that continue to use the selected PC files, and optionally residualized BEDs ending in `.residualized.bed.gz`.

## `workflows/common/calculate_phenotypePCs.wdl`

Workflow that computes phenotype PCs from any normalized molecular phenotype BED file.

**Steps:**
1. Runs `calculate_PCs.R` on the input BED file using the Gavish-Donoho method to select the number of PCs.

**Inputs:** Normalized BED file, output prefix, optional output suffix, resource parameters.

**Outputs:** Selected phenotype PCs TSV (`<OutputPrefix>_phenotype_PCs<OutputSuffix>.tsv`) and full rotated-PC TSV (`<OutputPrefix>_phenotype_PCs<OutputSuffix>.all.tsv`).

## `workflows/common/MergeCovariates.wdl`

Workflow that merges additional covariates, such as genotype PCs, and molecular phenotype PCs into a single covariate file ready for tensorQTL.

**Steps:**
1. Runs `MergeCovariates.R` to inner-join additional covariates and molecular PCs, then transpose the result.

**Inputs:** Additional covariates TSV with `sample_id`, Gavish-Donoho-selected molecular PCs TSV, output prefix, optional output suffix.

**Outputs:** Combined QTL covariate file (`<OutputPrefix>_QTL_covariates<OutputSuffix>.tsv`).

## `workflows/common/ResidualizePhenotypes.wdl`

One-task workflow component that residualizes a molecular phenotype BED and scales the residuals.

**Steps:**
1. Runs `ResidualizePhenotypes.R` on one normalized BED file.
2. If a covariate file is provided, residualizes each phenotype row against those covariates.
3. Centers and scales each residual row.

**Inputs:** Normalized BED file, optional merged covariates TSV, output BED filename, resource parameters.

**Outputs:** Residualized BED file, for example `<OutputPrefix>.expression.INT.residualized.bed.gz`.
