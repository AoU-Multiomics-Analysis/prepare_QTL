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

## MethBat chr21 pilot

Use [MethBat_jointseg_chromosome_subsets.wdl](MethBat_jointseg_chromosome_subsets.wdl)
with the existing legacy pb-CpG-tools manifest and indexes. The default settings
use chr21 and two disjoint donor groups, with nested subsets of 25, 50, and 100
donors per group. The default image uses a published digest from the
**MethBat Image** Actions run. It contains the PacBio v1.2.0 binary and the pilot
scripts. Use `methbat_docker` to override the default digest.

For each subset, the workflow converts the staged tracks to one unified 5mC
pileup per donor. The direct run uses the donor pileups. The upstream run uses
one cohort pileup and the same MethBat segmentation options.

- Keep a CpG if at least 90% of all selected donors have a combined call.
- Keep its haplotype signal if at least 60% of all selected donors have both
  haplotype calls. Order the two fractions from low to high within each donor
  before averaging the low and high values separately.
- Give each eligible donor equal weight at a site. Use the modified and
  unmodified counts from the input files. For a model-mode file, these are
  estimated counts. The upstream pileup stores averages on the 1/10,000
  fraction scale used by the v1.2.0 loader. This extra rounding can affect
  segmentation or output scores.

The workflow returns both region BEDs and both sets of score tracks, the cohort
pileup, CpG filter counts, and a paired method comparison for each run. The
comparison ignores file headers. It counts exact row matches, reports which
rows differ, and measures the maximum score difference on shared bedGraph
coordinates in percentage points. The existing boundary summaries still compare
the two disjoint groups at each size with the direct MethBat method.
