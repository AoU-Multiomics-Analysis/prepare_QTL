# Methylation workflow modules

These internal WDLs implement the public
[`AggregateMethylationCohort.wdl`](AggregateMethylationCohort.wdl) entry
point. They are grouped by pipeline stage so the public workflow contains only
cross-stage orchestration.

## Existing merged BED entry point

Use [`prepare_mQTL.wdl`](prepare_mQTL.wdl) when a cohort-wide methylation BED already exists. The first four input columns contain chromosome, start, end, and feature ID data. The remaining columns contain numeric methylation values, with one column per sample.

The sample list can be a headerless one-column file or use `sample_id`, `SampleID`, or `ID` as its header. Every requested sample must exist in the BED. The workflow keeps features with missingness strictly below `MissingnessThreshold` (5% by default), excludes random and non-autosomal contigs, and imputes retained missing values with each feature mean. At least two samples and two retained features are required for phenotype PC calculation. It writes raw, inverse-normalized, and scaled BED files. The raw BED keeps all selected samples. The INT and scaled BED files each use independent WGCNA connectivity filtering and have separate outlier reports. The workflow calculates PCs for both normalized files and can merge them with `AdditionalCovariates`.

| Module | Responsibility |
| --- | --- |
| `cohort_aggregation.wdl` | Expand the compact cohort manifest, validate cohort samples, merge each autosome, assemble cohort-wide pre-connectivity files and QC, and call the annotation workflow. |
| `connectivity.wdl` | Calculate preliminary phenotype PCs, build correlation covariates, analyze CpG correlation by chromosome, and remove connectivity outliers. |
| `annotation.wdl` | Own the cohort-filtered CpG annotation task and subworkflow called by `cohort_aggregation.wdl`. |
| `qtl_covariates.wdl` | Calculate final phenotype PCs and optionally merge additional covariates for TensorQTL. |

The root workflow owns the stable user-facing inputs and outputs. Add a task to
the module that owns its data, and expose only the files needed by another
stage.
