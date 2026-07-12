# Phenotype Normalization and Filtering

[Back to main README](../README.md)

This page summarizes the modality-specific normalization and filtering used by the molecular QTL prepare workflows. The prepare workflows write three BED matrices for each phenotype type:

- `.INT`: Rank-based inverse normal transformed values.
- `.scaled`: Centered and scaled values.
- `.raw`: Untransformed values after the modality-specific sample, feature, and BED-formatting filters.

Only `.INT` and `.scaled` outputs are used for phenotype PCs, optional covariate merging, optional residualization, and WGCNA connectivity outlier removal. Raw BEDs skip WGCNA outlier removal and keep all samples that remain after the initial sample-list filter.

## Shared WGCNA Connectivity Filtering

For `.INT` and `.scaled` matrices from expression, proteomics, and splicing, the prepare scripts remove sample connectivity outliers before writing the BED file:

1. Use the feature-by-sample matrix for the transformed phenotype branch.
2. Compute sample-sample biweight midcorrelations with `WGCNA::bicor(phenotype_matrix, use = "pairwise.complete.obs")`.
3. Convert correlations to adjacency with `0.5 + 0.5 * correlation`, replacing `NA` adjacency values with 0.
4. Calculate sample connectivity with `WGCNA::fundamentalNetworkConcepts()`.
5. Z-score sample connectivity values and remove samples with `Z_score < -3`.
6. Write removed samples to `<OutputPrefix>.<phenotype>.<normalization>.connectivity_outliers.tsv`.

The connectivity outlier TSV contains `SampleID` and `Z_score` columns. If a transformed matrix has fewer than 3 samples, fewer than 2 features, or zero/undefined connectivity variance, the scripts keep all samples and write an empty outlier TSV. Each outlier-removal round logs how many samples were removed and how many samples remain.

## Expression

Expression preparation is implemented in [`scripts/expression/PrepareExpression.R`](../scripts/expression/PrepareExpression.R) and wrapped by [`workflows/prepare_eQTL.wdl`](../workflows/prepare_eQTL.wdl).

**Input phenotype:** RNA-seq gene count GCT or TSV with `Name`, `Description`, and sample count columns.

**Filtering:**

- Keeps only sample columns listed in `--SampleList`.
- Filters genes before normalization, retaining genes with raw counts greater than 6 in at least 20% of retained samples.
- Extracts gene TSS positions from the GENCODE GTF using `type == "gene"` and strand-aware TSS coordinates.
- Inner-joins the expression matrix to the GTF-derived TSS table by `gene_id`; genes without matching TSS annotation are dropped from the BED.
- Removes WGCNA connectivity outlier samples from `.INT` and `.scaled` matrices only.

**Normalization and transforms:**

- Runs edgeR TMM normalization with `edgeR::calcNormFactors()`.
- Converts TMM-normalized counts to CPM with `edgeR::cpm(..., log = FALSE)`.
- `.raw`: TMM-normalized CPM values, without rank-normalization, log2 transform, centering, scaling, or WGCNA outlier removal.
- `.INT`: applies `RNOmni::RankNorm()` to each gene across samples using unlogged CPM values.
- `.scaled`: applies `log2(CPM + 1)` to each gene across samples, then centers and scales with `scale(..., center = TRUE, scale = TRUE)`.

## Proteomics

Proteomics has an optional Olink preprocessing workflow in [`scripts/proteomics/NormalizeProteomics.R`](../scripts/proteomics/NormalizeProteomics.R), followed by pQTL BED preparation in [`scripts/proteomics/PrepareProteomics.R`](../scripts/proteomics/PrepareProteomics.R) and [`workflows/prepare_pQTL.wdl`](../workflows/prepare_pQTL.wdl).

### Olink NPX Preprocessing

[`NormalizeProteomics.R`](../scripts/proteomics/NormalizeProteomics.R) reads Olink NPX parquet files and writes a filtered long-format table that can be passed to `PrepareProteomics.R`.

**Filtering before normalization:**

- Reads all parquet files in `--OlinkDataDir` with `OlinkAnalyze::read_NPX()`.
- Requires `SampleID`, `PlateID`, `AssayType`, `OlinkID`, `NPX`, and `UniProt`; downstream QC filtering also uses `Assay`, `AssayQC`, and `SampleQC`.
- Keeps assay rows and removes control samples matching `CONTROL_SAMPLE`, `Control`, or `CONT`.
- Requires the requested `--ReferencePlate` to be present.

**Normalization:**

- Calculates reference medians per `OlinkID` from the reference plate.
- Runs `OlinkAnalyze::olink_normalization()` separately by plate using those reference medians.
- Writes the full median-normalized table to `<OutputPrefix>_median_normalized.tsv.gz`.

**Filtering after normalization:**

- Removes assays with `AssayQC != "PASS"`.
- Removes samples with `SampleQC != "PASS"` when the warning/failure is not already explained by a removed assay.
- Removes control samples and keeps assay rows.
- Strips the plate suffix from `SampleID` with `str_remove(SampleID, "_.*")`.
- For samples appearing on multiple plates, keeps the first observed sample-plate pair; samples appearing on one plate keep that pair.
- Removes any sample with missing `NPX`.
- Removes UniProt IDs `P32455` and `Q02750`.
- Writes `<OutputPrefix>_npx_values.tsv.gz` with `ResearchID`, `UniProt`, and `PCNormalizedNPX`.

### pQTL BED Preparation

**Input phenotype:** normalized long-format proteomics data with `PCNormalizedNPX`, `UniProt`, and either `SampleID` or `ResearchID`.

**Filtering:**

- Selects `SampleID` or `ResearchID` as the sample identifier, choosing the column with the largest overlap with `--SampleList`.
- Keeps only rows whose selected sample ID appears in `--SampleList`.
- Maps UniProt IDs to Ensembl gene IDs with BioMart.
- Extracts gene TSS positions from the GENCODE GTF.
- Joins UniProt-to-gene mappings to TSS locations by `gene_id`.
- Drops proteins without TSS annotation after the join.
- Removes WGCNA connectivity outlier samples from `.INT` and `.scaled` matrices only.

**Normalization and transforms:**

- Uses `PCNormalizedNPX` as the starting value for pQTL BED preparation.
- `.raw`: `PCNormalizedNPX` values, without rank-normalization, centering, scaling, or WGCNA outlier removal.
- `.INT`: applies `RNOmni::RankNorm()` to each protein across samples.
- `.scaled`: centers and scales each protein across samples with `scale(..., center = TRUE, scale = TRUE)`.

## Splicing

Splicing preparation is implemented in [`scripts/splicing/PrepareSpliceData.R`](../scripts/splicing/PrepareSpliceData.R) and wrapped by [`workflows/prepare_sQTL.wdl`](../workflows/prepare_sQTL.wdl).

**Input phenotype:** LeafCutter-style BED with the first four columns containing interval and phenotype identifiers, followed by sample columns.

**Filtering:**

- Keeps the first four BED columns and sample columns listed in `--SampleList`.
- Does not apply an additional feature-level expression or missingness filter in the prepare script.
- Extracts an `ENSG` ID from the phenotype ID for sorting.
- Sorts output rows by chromosome, start, end, and extracted Ensembl gene ID.
- Removes WGCNA connectivity outlier samples from `.INT` and `.scaled` matrices only.

**Normalization and transforms:**

- Uses the input splice values as the starting matrix.
- `.raw`: input splice values after sample selection and BED formatting, without rank-normalization, centering, scaling, or WGCNA outlier removal.
- `.INT`: applies `RNOmni::RankNorm()` to each splice phenotype across samples.
- `.scaled`: centers and scales each splice phenotype across samples with `scale(..., center = TRUE, scale = TRUE)`.
