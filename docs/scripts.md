# R Scripts

[Back to main README](../README.md)

All scripts are written in R and are invoked from the command line with `Rscript`. They are bundled in the Docker image defined in [`envs/PhenotypePCs/Dockerfile`](../envs/PhenotypePCs/Dockerfile). For modality-specific normalization and filtering details, see [Phenotype normalization and filtering](phenotype-normalization-filtering.md).

## Shared Dual-Output Behavior

[`PrepareExpression.R`](../scripts/expression/PrepareExpression.R), [`PrepareProteomics.R`](../scripts/proteomics/PrepareProteomics.R), and [`PrepareSpliceData.R`](../scripts/splicing/PrepareSpliceData.R) each write three BED files:

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

## `scripts/expression/PrepareExpression.R`

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

## `scripts/proteomics/PrepareProteomics.R`

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

## `scripts/proteomics/NormalizeProteomics.R`

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
- `<OutputPrefix>_npx_values.tsv.gz`: Filtered long-format table with `ResearchID`, `UniProt`, and `PCNormalizedNPX`. This file can be used as `--ProteomicData` for `PrepareProteomics.R` or as `ProteomicData` for [`workflows/proteomics/prepare_pQTL.wdl`](../workflows/proteomics/prepare_pQTL.wdl).

## `scripts/splicing/PrepareSpliceData.R`

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

## Methylation scripts

The pb-CpG-tools 5mC workflow is split into shard filtering, per-autosome cohort merging, and final aggregation stages. These scripts source `scripts/methylation/MethylationUtils.R` for BED parsing, validation, QC metrics, transformations, and plotting helpers.

- `scripts/methylation/FilterMethylationShard.R` applies chromosome, minimum-coverage, and extreme-coverage QC to a manifest shard. It writes one QC-flagged call table per autosome plus one-row-per-sample QC summaries; it does not create a redundant filtered-call copy.
- `rust/methylation_filter` provides the equivalent streaming parser used only by `ProcessMethylationSample.wdl`. It makes two sequential passes over the localized per-sample BED (to determine the sample-specific coverage fence, then emit chromosome files), so it avoids materializing the full call table in memory.
- `rust/methylation_merge` is the bounded-memory, per-chromosome cohort merger used by `AggregateMethylationCohort.wdl`. It k-way merges coordinate-sorted per-sample call files, aggregates one CpG at a time, and creates staged compressed runs when the cohort has more than 128 samples.
- `scripts/methylation/BuildMethylationCohortSamples.R` validates sample IDs, constructs the cohort sample order, and combines sample-QC files once before the per-chromosome cohort scatter.
- `scripts/methylation/AnalyzeMethylationCpGCorrelation.R` streams a sorted chromosome INT methylation BED and computes covariate-adjusted local Pearson correlations. It writes cluster summaries, QC plots, and one most-connected representative per correlated cluster (plus singleton representatives); it does not remove CpGs from the phenotype BED.
- `scripts/methylation/BuildMethylationCorrelationCovariates.R` formats preliminary phenotype PCs and optional additional covariates in TensorQTL orientation for the correlation stage. Additional-covariate rows outside the phenotype-PC cohort are discarded; all phenotype-PC samples must be represented.
- `scripts/methylation/FinalizeMethylationConnectivity.R` uses a capped representative-CpG set for landmark sample connectivity and removes failing samples consistently from final raw/INT BEDs and filtered calls.
- `scripts/methylation/MergeMethylationCohort.R` reduces one chromosome's shard outputs, derives passing calls from `per_sample_qc_pass`, applies cohort sample-presence and MAD filters, records a sample-normalized methylation–coverage correlation diagnostic, mean-imputes retained features, and writes chromosome-level raw and INT QTL phenotype BEDs.
- `scripts/methylation/AggregateMethylationChromosomes.R` writes the final sample-QC table and creates global filter summaries and plots after the WDL has concatenated chromosome-level tables as header-aware compressed streams.
- `scripts/methylation/AnnotateMethylationSites.R` annotates retained sites with the nearest strand-aware TSS, promoter/gene-body/intergenic and exon/intron/CDS/UTR context, overlapping ENCODE cCREs, and UCSC CpG-island, shore, shelf, or open-sea context.

See the [PacBio 5mC QTL workflow guide](methylation-qtl.md) for the input schema, all command-line options, QC logic, outputs, and QTL phenotype format.

## `scripts/common/calculate_PCs.R`

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

## `scripts/common/MergeCovariates.R`

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

## `scripts/common/ResidualizePhenotypes.R`

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
