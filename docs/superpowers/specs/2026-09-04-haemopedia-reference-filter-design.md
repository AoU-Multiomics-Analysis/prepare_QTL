# Haemopedia comparison and cell-type-specific expression filtering

## Objective and scope

Add post-export filtering to the maintained cell-type expression workflow in
`prepare_QTL`. Produce nonnegative, reference-supported phenotype BEDs for eQTL
preparation, with a complete record of exclusions and reference comparisons.
Keep original TCA estimates and the existing numerically cleaned TCA model.
Do not refit TCA, change its statistical model, clip estimates, or transform the
values stored in the filtered BEDs.

The user approved removing a gene only from the affected cell type, the reference
mapping below, and equal weighting of reference samples. This document specifies
the remaining interface and edge-case rules for review before implementation.

## Starting evidence

- The workflow exports TCA estimates on the linear CPM scale.
- `PrepareExpression.R` rejects negative values supplied through `CpmBed`.
- The existing summary computes statistics on stored CPM-scale values. It does
  not compute the mean of sample-level log expression.
- The user supplied Haemopedia raw-count processing: `edgeR::DGEList`,
  `edgeR::calcNormFactors`, then `edgeR::cpm(log = FALSE)`.
- The supplied reference header has 42 expression columns covering 12 sorted
  populations. The raw file itself has not been inspected locally.
- These are human data. Reference donor matching to the deconvolution cohort is
  not established; the comparison is across genes, not paired donors.

Reference: [Haemopedia sample descriptions](https://www.haemosphere.org/datasets/samples?datasetName=Haemopedia-Human-RNASeq).

## Workflow placement and inputs

Keep the existing fitting, numerical model cleanup, export, raw gene summary,
and reconstruction QC steps. Add reference preparation and a post-export filter
stage. Feed the filtered BEDs and their inventory into the existing eQTL scatter.
Keep all stages within the maintained cell-type workflow directory.

Proposed public inputs:

- `File? haemopedia_counts`: the raw Haemopedia count table, plain or gzip TSV.
  When omitted, run negative-value filtering only and record
  `reference_not_provided`. This preserves use without an external reference.
- `Float reference_min_mean_log2_cpm1 = 0.01`: strict lower expression threshold
  for both datasets when a reference is supplied.
- `Float? reference_residual_cutoff`: positive cutoff for absolute internally
  standardized regression residuals. Omitted means report residuals without
  removing genes based on residuals. No data-derived automatic cutoff.

Specifying a residual cutoff without a reference is an input error. These inputs
are separate from the HSPE pseudocount, which is not used in this comparison.

## Reference preparation

Read the first column as gene IDs and all remaining columns as samples. Accept
the source's empty first-column header or `V1`, and an explicit `gene_id` header;
do not require a hand-edited source file. Validate unique nonempty IDs and sample
names, finite nonnegative counts, and positive library totals. Do not treat a
normalized CPM or TPM matrix as raw counts under this input contract.

Apply TMM normalization and linear CPM calculation to the full reference count
matrix before restricting its gene rows to the deconvolution gene set. Do not
normalize each cell type separately or recalculate CPM after filtering. Record
library sizes, normalization factors, package versions, and resolved sample map.

Remove only the terminal replicate suffix (`\\.[0-9]+$`) from sample names, then
use the explicit mapping below. Unknown sample prefixes are errors rather than
guessed cell types. Require all mapped subpopulations for each compared lineage;
do not silently build a partial reference when an expected subtype is absent.

| Pipeline cell type | Reference populations |
| --- | --- |
| B cells | NveB, MemB |
| CD4 T cells | CD4T |
| CD8 T cells | CD8T |
| NK cells | NK |
| Monocyte/myeloid | Mono, MonoNonClassical |
| Neutrophils | Neut |
| Eosinophils | Eo |
| Dendritic cells | myDC, myDC123, pDC |

Pool sample columns within each mapped lineage and give each sample equal
weight. Do not average subtype means with equal subtype weights. Thus five
naive-B samples contribute more than three memory-B samples. Reference donors
may contribute samples from multiple subtypes; this is a descriptive reference
profile, not an independent-replicate statistical test.

Mast cells and gamma-delta T cells have no reference match. Apply their negative
filter only and record `no_reference_cell_type`; do not assign a substitute.

Match genes by Ensembl gene ID. Allow removal of an anchored numeric version
suffix from Ensembl gene IDs for matching only. Preserve the original BED IDs.
Fail on duplicate matching keys after normalization; do not silently aggregate
them. Gene symbols are not an implicit fallback. Record unmatched gene counts.

## Per-cell-type filter and report sequence

Process BEDs in bounded row chunks. Do not load the complete cohort tensor into
memory or convert it to a full sample-by-gene long table.

1. Validate BED coordinates, unique gene IDs, sample names and order, and finite
   numeric values. Missing or nonfinite values are input errors, not silently
   omitted measurements.
2. Before filtering, record each gene's negative count, negative percentage,
   minimum CPM, mean negative CPM, and number of samples. Mean negative CPM is
   missing when the count is zero. Negative means strictly less than zero, with
   no tolerance or clipping. Exclude a gene from this BED if any sample is
   negative. Do not exclude the gene from other cell types for this reason.
3. For eligible nonnegative genes, calculate `mean_log2_cpm1` and
   `median_log2_cpm1` from sample-level `log2(CPM + 1)`. This is not
   `log2(mean(CPM) + 1)`. Calculate the same statistics in the reference. Keep
   separate CPM-scale summaries where useful; never label log statistics CPM.
4. For a mapped reference lineage, exclude unmatched genes and require mean log
   expression strictly greater than the threshold in both datasets. A missing
   reference gene is not a measured zero. Keep separate reason flags for
   negative values, reference absence, and low expression in either dataset.
5. On the remaining shared genes, fit an ordinary least-squares regression with
   intercept: `deconvolution_mean_log2_cpm1 ~ reference_mean_log2_cpm1`.
6. Record Pearson r, Spearman rho, regression R-squared, slope, intercept, sample
   counts for both profiles, and number of genes. Save fitted values, signed
   residuals, and internally standardized residuals for every eligible gene.
   The standardized residual is `e / (s * sqrt(1 - h))`, using residual standard
   error `s` and leverage `h` from this single baseline fit.
7. If a cutoff was supplied, remove genes with absolute standardized residual
   strictly greater than that cutoff. Use one pass; do not repeatedly refit and
   remove more genes. Report metrics on retained genes separately, with a new
   descriptive fit, without replacing the original metrics or residuals.
8. Write retained rows with unchanged linear CPM values, sample order, gene IDs,
   and BED coordinates. Different cell types may retain different genes. Keep
   an exclusion table for all original gene-cell-type pairs, including columns
   that identify the final retention decision and all evaluated reason flags.

At least three shared genes and variation in both mean-expression vectors are
required for correlation/regression reporting. Report unavailable metrics with
an explicit reason, not zero. Undefined standardized residuals, including a
zero residual scale or unit leverage, are recorded as unavailable. If residual
filtering is requested and those residuals cannot be calculated, fail with a
clear diagnostic rather than silently skipping the requested filter.

If a cell type has no genes after mandatory filtering, stop before eQTL
preparation with its cell type and filter counts in the log. Do not send an empty
BED downstream. Existing eQTL checks still apply, including constant-gene checks.

## Outputs and plots

Preserve existing original BED and original gene-summary outputs. Add explicitly
named filtered BEDs and an inventory with per-cell-type retained gene counts.
Do not enforce identical gene coordinates across different filtered BEDs: the
existing raw summary's cross-cell-type equality check remains on original BEDs.

Add:

- A gene-by-cell-type negative-estimate summary and exclusion/comparison table.
- A reference mean-log-expression table and sample-mapping/normalization report.
- Per-cell-type before/after filter counts and comparison metrics.
- A percentage-of-genes-with-negatives plot and negative-fraction heatmap.
  For a large heatmap, show at most 100 genes ranked by maximum negative fraction
  across cell types, with gene ID as the tie-breaker. Keep every gene in the table.
- Per-cell-type scatter plots with reference on x and deconvolution on y, an
  identity line and regression line, and residual-exclusion status shown when
  applicable. Include a residual plot. Do not silently downsample genes.

Use clear axis labels, a minimal theme, and no plot titles or subtitles. Keep
cell types identifiable through file names, facets, or axis labels. Use separate
plot files for the pre-filter and post-residual-filter comparisons.

Extend the final QTL manifest without removing existing columns. Keep its Terra
entity/cell-type IDs and existing QTL outputs. Add source and filtered CPM BED
cloud paths and links to the shared comparison/filter reports. These links may
repeat across rows when a report contains all cell types. Keep File localization
separate from String cloud-path reporting as in the current manifest fix.

## Scientific limits and claims

The filter identifies genes that meet a specified nonnegativity and reference
expression rule; it does not prove accurate cell-specific recovery. The B-cell
model group includes plasma cells, absent from this reference. The myeloid model
group includes macrophage components, while the reference contains monocytes.
The dendritic comparison combines distinct sorted populations. Record these
limitations with the mapping.

Residual differences can reflect subtype composition, biological differences,
annotation or measurement differences, and deconvolution error. Correlation
across genes is descriptive and is not donor-level validation. For regression
with one predictor and an intercept, R-squared equals Pearson r squared; high
correlation does not establish agreement with the identity line. Reference-based
selection changes the tested gene universe and can remove genuine biology.
Any improved agreement after residual filtering is selection-dependent, not an
independent validation result.

## Implementation and validation requirements

Use R/tidyverse patterns, existing task logging, WDL 1.0, and existing global
retry settings. Create serialized files only in task scope, never at workflow
scope. Use the pinned container dependencies and GitHub Actions for container
builds and smoke tests; do not build a local Docker image for this change.

Tests must cover:

- Exact TMM CPM parity with the supplied edgeR sequence; normalize before subset.
- Mean of logs versus log of mean; unequal subtype sample counts and suffix map.
- A gene negative in one cell type but retained in another; zeros and tiny
  negative values; preservation of every sample and coordinate in retained rows.
- Strict expression threshold equality, missing reference genes, unsupported
  cell types, omitted reference, and duplicate/versioned identifier handling.
- Known regression/correlation results, fixed one-pass residual exclusions,
  perfect-fit and too-few-gene edge cases, and unchanged stored CPM values.
- Filtered BED routing into eQTL, unequal gene sets across BEDs, manifest cloud
  paths, and retention of all previous public outputs.
- Static workflow-scope file-writing regression, WDL validation, and a small
  GitHub Actions end-to-end fixture through filtered BED and QTL manifest output.

No cloud jobs are authorized by this design. Report local and GitHub checks
separately from a complete Terra run; the new workflow is not yet tested on Terra.
