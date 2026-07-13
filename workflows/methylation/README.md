# Methylation workflow modules

These internal WDLs implement the public
[`AggregateMethylationCohort.wdl`](AggregateMethylationCohort.wdl) entry
point. They are grouped by pipeline stage so the public workflow contains only
cross-stage orchestration.

| Module | Responsibility |
| --- | --- |
| `cohort_aggregation.wdl` | Expand the compact cohort manifest, validate cohort samples, merge each autosome, assemble cohort-wide pre-connectivity files and QC, and call the annotation workflow. |
| `connectivity.wdl` | Calculate preliminary phenotype PCs, build correlation covariates, analyze CpG correlation by chromosome, and remove connectivity outliers. |
| `annotation.wdl` | Own the cohort-filtered CpG annotation task and subworkflow called by `cohort_aggregation.wdl`. |
| `qtl_covariates.wdl` | Calculate final phenotype PCs and optionally merge additional covariates for TensorQTL. |

The root workflow owns the stable user-facing inputs and outputs. Add a task to
the module that owns its data, and expose only the files needed by another
stage.
