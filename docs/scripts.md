# R Scripts

[Back to main README](../README.md)

All scripts are written in R and are invoked from the command line with `Rscript`. They are bundled in the Docker image defined in [`envs/PhenotypePCs/Dockerfile`](../envs/PhenotypePCs/Dockerfile). For modality-specific normalization and filtering details, see [Phenotype normalization and filtering](phenotype-normalization-filtering.md).

## Shared Dual-Output Behavior

[`PrepareExpression.R`](../scripts/PrepareExpression.R), [`PrepareProteomics.R`](../scripts/PrepareProteomics.R), and [`PrepareSpliceData.R`](../scripts/PrepareSpliceData.R) each write three BED files:

- `.INT`: Rank-based inverse normal transformed molecular phenotypes.
- `.scaled`: Molecular phenotypes transformed with `scale(..., center = TRUE, scale = TRUE)`. Expression CPMs are transformed with `log2(CPM + 1)` before centering/scaling; proteomics and splicing values are centered/scaled directly.
- `.raw`: Untransformed phenotype values after sample/feature filtering and BED formatting.

For each `.INT` and `.scaled` matrix, the scripts compute WGCNA sample connectivity outliers from the feature-by-sample matrix before writing the BED. Raw BEDs keep all samples after the initial sample-list filter.

WGCNA outlier detection uses the same implementation in all three prepare scripts:

1. Compute sample-sample biweight midcorrelations with `WGCNA::bicor(phenotype_matrix, use = "pairwise.complete.obs")`.
2. Convert correlations to adjacency with `0.5 + 0.5 * correlation`, replacing `NA` adjacency values with 0.
3. Calculate sample connectivity with `WGCNA::fundamentalNetworkConcepts()`.
4. Z-score the connectivity values and remove samples with `Z_score < -3`.
5. Write removed samples to `<OutputPrefix>.<modality>.<normalization>.connectivity_outliers.tsv` with `SampleID` and `Z_score` columns.

If a transformed matrix has fewer than 3 samples, fewer than 2 features, or zero/undefined connectivity variance, the scripts keep all samples and write an empty outlier TSV. Each outlier-removal round logs how many samples were removed and how many remain; raw outputs log that outlier removal was skipped.

The scripts still accept `--RankNormalize` for backwards compatibility, but the option is deprecated and no longer selects a single output mode.

## `scripts/PrepareExpression.R`

Prepares RNA-seq gene expression data for eQTL analysis.

**What it does:**
1. Loads raw count data from a GCT-formatted file (`.gct` or tab-separated).
2. Loads a GENCODE GTF annotation file and extracts transcription start site (TSS) positions for each gene.
3. Filters to a specified list of samples.
4. Filters genes by expression count, retaining genes with counts > 6 in at least 20% of samples.
5. Normalizes counts using edgeR TMM normalization followed by CPM transformation.
6. Creates rank-based inverse normal transformed, `log2(CPM + 1)` centered/scaled, and raw CPM matrices. The log2 transform is only applied to the centered/scaled expression branch, not the INT branch.
7. Merges each matrix with TSS locations and writes compressed BED files.

**Inputs:**
- `--CountGCT`: GCT or TSV file of raw RNA-seq count data.
- `--AnnotationGTF`: GENCODE GTF file used to extract TSS locations.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Deprecated compatibility option. Both `.INT` and `.scaled` outputs are always written.

**Outputs:**
- `<OutputPrefix>.expression.INT.bed.gz`
- `<OutputPrefix>.expression.scaled.bed.gz`
- `<OutputPrefix>.expression.raw.bed.gz`
- `<OutputPrefix>.expression.INT.connectivity_outliers.tsv`
- `<OutputPrefix>.expression.scaled.connectivity_outliers.tsv`

## `scripts/PrepareProteomics.R`

Prepares Olink proteomics data for pQTL analysis.

**What it does:**
1. Reads Olink proteomics data from a TSV, gzipped TSV, or Parquet file.
2. Loads a GENCODE GTF file to extract TSS positions.
3. Uses Ensembl BioMart to map UniProt IDs to Ensembl gene IDs.
4. Filters to a specified list of samples using either `SampleID` or `ResearchID`.
5. Pivots data to wide format with proteins as columns and samples as rows.
6. Creates rank-based inverse normal transformed, centered/scaled, and raw protein matrices.
7. Joins the UniProt-to-Ensembl mapping to TSS locations by `gene_id` and writes compressed BED files.

**Inputs:**
- `--ProteomicData`: TSV, TSV.gz, or Parquet file of normalized Olink protein expression data. The file must contain `PCNormalizedNPX`, `UniProt`, and at least one sample ID column: `SampleID` or `ResearchID`. If both sample ID columns are present, the script uses the one with the largest overlap with `--SampleList`.
- `--AnnotationGTF`: GENCODE GTF file used to extract TSS locations.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Deprecated compatibility option. Both `.INT` and `.scaled` outputs are always written.

**Outputs:**
- `<OutputPrefix>.protein.INT.bed.gz`
- `<OutputPrefix>.protein.scaled.bed.gz`
- `<OutputPrefix>.protein.raw.bed.gz`
- `<OutputPrefix>.protein.INT.connectivity_outliers.tsv`
- `<OutputPrefix>.protein.scaled.connectivity_outliers.tsv`

## `scripts/NormalizeProteomics.R`

Median-normalizes Olink NPX parquet files and writes filtered long-format protein values for pQTL preparation.

**What it does:**
1. Loads all parquet files in `--OlinkDataDir` with `OlinkAnalyze::read_NPX`.
2. Keeps Olink assay rows and removes control samples.
3. Calculates assay-level reference medians from `--ReferencePlate`.
4. Runs Olink plate normalization with `OlinkAnalyze::olink_normalization`.
5. Writes the full median-normalized Olink table.
6. Removes assays with non-`PASS` assay QC, removes samples with non-`PASS` sample QC when the warning/failure is not already explained by a removed assay, filters missing samples, keeps one valid sample-plate row per sample, removes UniProt IDs `P32455` and `Q02750`, and writes the long-format pQTL input table.

**Inputs:**
- `--OlinkDataDir`: Directory containing Olink NPX parquet files. Input files must include `SampleID`, `PlateID`, `AssayType`, `OlinkID`, `NPX`, and `UniProt`; QC filtering also uses `Assay`, `AssayQC`, and `SampleQC`.
- `--OutputPrefix`: Prefix for output files.
- `--OutputDir`: Directory for output files.
- `--ReferencePlate`: Plate ID used to calculate reference medians.

**Outputs:**
- `<OutputPrefix>_median_normalized.tsv.gz`: Full median-normalized Olink table.
- `<OutputPrefix>_npx_values.tsv.gz`: Filtered long-format table with `ResearchID`, `UniProt`, and `PCNormalizedNPX`. This file can be used as `--ProteomicData` for `PrepareProteomics.R` or as `ProteomicData` for [`workflows/prepare_pQTL.wdl`](../workflows/prepare_pQTL.wdl).

## `scripts/PrepareSpliceData.R`

Prepares LeafCutter splice junction data for sQTL analysis.

**What it does:**
1. Loads a BED file of splice junction quantifications produced by LeafCutter.
2. Filters columns to a specified list of samples.
3. Creates rank-based inverse normal transformed, centered/scaled, and raw splice matrices.
4. Sorts junctions by chromosome, start, end, and Ensembl gene ID.
5. Writes the normalized splice data to compressed BED files.

**Inputs:**
- `--SpliceData`: BED file of splice junction quantifications from LeafCutter.
- `--SampleList`: File containing the list of sample IDs to include.
- `--OutputPrefix`: Prefix for the output file.
- `--RankNormalize`: Deprecated compatibility option. Both `.INT` and `.scaled` outputs are always written.

**Outputs:**
- `<OutputPrefix>.splicing.INT.bed.gz`
- `<OutputPrefix>.splicing.scaled.bed.gz`
- `<OutputPrefix>.splicing.raw.bed.gz`
- `<OutputPrefix>.splicing.INT.connectivity_outliers.tsv`
- `<OutputPrefix>.splicing.scaled.connectivity_outliers.tsv`

## `scripts/MergeMethylationCalls.R`

Merges per-sample pb-CpG-tools 5mC calls after per-sample and cohort-level coverage QC. Point each manifest `file_path` at the sample's `.combined.bed.gz` output. The script safely skips pb-CpG-tools metadata lines and requires its `#chrom`, `begin`, `end`, `mod_score`, `type`, and `cov` columns.

**What it does:**
1. Removes contigs matching `--FilterChroms` (by default `X|Y|M|_`).
2. Calculates each sample's median coverage and uses a Tukey far-out fence on `log10(cov)` to flag extreme coverage calls.
3. Removes calls below `--MinCoverage` and all extreme-coverage calls.
4. Retains only sites that pass per-sample QC in at least `max(ceiling(n_samples * MinSampleFraction), MinSamples)` samples.
5. Writes the retained long-format calls and per-sample/per-site QC tables. It also writes all-site metadata, including coverage and methylation mean, standard deviation, coefficient of variation (CV), coverage fractions, and the final `keep_site` flag. When `--ValueColumn` is supplied, it is used for methylation metrics and for the site-by-sample matrix.

The run log reports, for every sample, input sites, chromosome-filtered sites, sites failing `MinCoverage`, sites failing the extreme-coverage filter after passing `MinCoverage`, and sites passing both. It also reports cohort-wide site counts after each stage and the numbers passing/failing the final sample-presence threshold.

For a sharded run, use `--PerSampleOnly` for each shard. It writes both `<OutputPrefix>.methylation.per_sample_filtered.long.tsv.gz` and `<OutputPrefix>.methylation.per_sample_qc.long.tsv.gz`, plus sample QC, without applying the cohort-level threshold. Run the final aggregation with `--FilteredCallList`, `--AllCallList`, `--TotalSamples`, and optionally `--FilteredSampleQcList`; only this final aggregation applies `--MinSampleFraction` and `--MinSamples` across the complete cohort.

**Inputs:**
- `--InputManifest`: TSV with `sample_id` and `file_path` columns.
- `--FilteredCallList`: In final sharded-merge mode, a one-path-per-line list of per-sample-QC-passing call files.
- `--AllCallList`: In final sharded-merge mode, a one-path-per-line list of all chromosome-filtered calls annotated with per-sample QC flags; required to calculate all-site metadata.
- `--OutputPrefix`: Prefix for output files.
- `--MinCoverage`: Minimum coverage required for a call (default `10`).
- `--MinSampleFraction`: Minimum cohort fraction required for a site (default `0.8`).
- `--MinSamples`: Optional additional minimum count of passing samples (default `0`).
- `--FilterChroms`: Regex for contigs to remove; pass an empty string to retain all contigs (default `X|Y|M|_`).
- `--FenceK`: Tukey far-out-fence multiplier (default `3`).
- `--ValueColumn`: pb-CpG-tools methylation column used for all-site methylation metrics and the QTL BED (default `mod_score`).
- `--ValueMultiplier`: Multiplier for `ValueColumn` before QTL output (default `0.01`, converting pb-CpG `mod_score` percent to 0–1 beta values).

**Outputs:**
- `<OutputPrefix>.methylation.filtered.long.tsv.gz`
- `<OutputPrefix>.methylation.site_qc.tsv.gz`
- `<OutputPrefix>.methylation.site_metadata.tsv.gz`: One row for every site observed after chromosome filtering, including raw-call and per-sample-QC-passing coverage/methylation metrics and `keep_site`.
- `<OutputPrefix>.methylation.sample_qc.tsv`
- `<OutputPrefix>.methylation.matrix.bed.gz`: TensorQTL-compatible beta-value BED, with `#chr`, `start`, `end`, and `phenotype_id` as its first four columns.

## `scripts/calculate_PCs.R`

Computes phenotype principal components (PCs) from a normalized BED file.

**What it does:**
1. Reads a normalized BED file, such as expression, splicing, or proteomics.
2. Runs PCA on the molecular phenotype data using `PCAtools::pca`.
3. Selects the optimal number of PCs using the Gavish-Donoho method (`chooseGavishDonoho`).
4. Writes the phenotype PCs to a TSV file.

**Inputs:**
- `--bed_file`: A normalized BED file of molecular phenotype data.
- `--output_prefix`: Prefix for the output file.
- `--output_suffix`: Optional suffix added before the `.tsv` extension.

**Output:** `<output_prefix>_phenotype_PCs<output_suffix>.tsv`

## `scripts/MergeCovariates.R`

Merges additional covariates, such as genotype PCs, and molecular phenotype PCs into a single covariate file for QTL analysis.

**What it does:**
1. Reads additional covariates and molecular phenotype PCs from TSV files.
2. Inner-joins the two tables by sample ID, using `sample_id` in the additional covariates file and `ID` in the phenotype PCs file.
3. Transposes the merged table so that rows are covariate names and columns are sample IDs, as required by tensorQTL.
4. Writes the combined covariate matrix to a TSV file.

**Inputs:**
- `--GenotypePCs`: TSV file of additional covariates. This file must contain a `sample_id` column.
- `--MolecularPCs`: TSV file of molecular phenotype PCs. This file must contain an `ID` column.
- `--OutputPrefix`: Prefix for the output file.
- `--OutputSuffix`: Optional suffix added before the `.tsv` extension.

**Output:** `<OutputPrefix>_QTL_covariates<OutputSuffix>.tsv`

## `scripts/ResidualizePhenotypes.R`

Residualizes a normalized molecular phenotype BED against merged covariates and then centers/scales each row of residuals.

**What it does:**
1. Reads a normalized molecular phenotype BED file.
2. Optionally reads merged tensorQTL-style covariates with covariates as rows and samples as columns.
3. Aligns samples between the BED and covariate matrix.
4. For each phenotype row, regresses phenotype values on the covariate design matrix.
5. Centers and scales each residual row across samples.
6. Writes a BED file with the original first 4 metadata columns and residualized/scaled sample values.

If no covariates are supplied, the script skips residualization and only centers/scales the input phenotype rows.

**Inputs:**
- `--InputBed`: Normalized molecular phenotype BED file.
- `--Covariates`: Optional merged covariates TSV in tensorQTL format.
- `--OutputFile`: Optional output BED file. If omitted, the script inserts `.residualized` before `.bed.gz`.

**Output:** Residualized and centered/scaled BED file, for example `<OutputPrefix>.expression.INT.residualized.bed.gz`.
