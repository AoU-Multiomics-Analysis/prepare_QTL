# Cell-type-specific expression workflows

These WDL 1.0 workflows prepare cell-type-specific expression from human
whole blood. `prepare_QTL` is the only maintained and canonical source for
this implementation.

## Public workflows

Use
[`deconvolution.wdl`](../workflows/cell_type_specific_expression/deconvolution.wdl)
to run hspe, process the proportions, fit TCA, and export one BED for each
retained major cell type.

Use
[`prepare_cell_type_eQTL.wdl`](../workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl)
to run the same deconvolution and then prepare QTL inputs for each retained
cell type. This workflow calls the existing single-matrix
[`prepare_eQTL.wdl`](../workflows/expression/prepare_eQTL.wdl). It does not
copy the expression-QTL tasks.

## HSPE installation and migration

The image installs HSPE 0.1 from the
[upstream source archive](https://gjhunt.github.io/hspe/hspe_0.1.tar.gz), not CRAN.
It verifies SHA-256
`f6c7c27ba2aec9b042cb8f793d7b7565023461812cdebdc0228b240d67b68c78`
before installation. The micromamba environment supplies pinned dependencies:
DEoptimR 1.1-4 and nloptr 2.2.1.

HSPE replaces dtangle for proportion estimation. Rename the `dtangle_*` inputs
to `hspe_*` in existing Terra configurations or input JSON files. This includes
`marker_fraction`, `quantile_normalize`, `cpu`, `memory`, and `disk_gb` suffixes.
The estimator outputs use the `hspe` prefix. The
`estimated_proportions` output name and the precomputed-proportion schema do not
change. Use the matching workflow and container release; an older container
does not contain the batch preparation, fitting, and merge scripts.

The estimator uses HSPE's default DEoptimR optimizer and stopping criteria.
`random_seed` controls HSPE and TCA independently. HSPE 0.1 fits samples
sequentially within each batch; Terra can run separate batches concurrently.
The metadata records the HSPE version, optimizer version, seed, markers, and
LM22 overlap. Marker selection remains `ratio` with fraction `0.10` by default.
No new gene-expression filter or normalization is introduced by this migration.

## HSPE sample batches

`PrepareHspeBatches` reads the gene-type-filtered BED and GTF once for HSPE.
It maps gene IDs, combines duplicate symbols in linear CPM, applies the log2
transform, and intersects with LM22. If joint quantile normalization is enabled,
it runs once on all reference and bulk profiles before marker selection.
The task selects ratio markers once from the reference and writes small RDS
files that contain only selected-marker expression for each sample batch.
Workers do not receive the full BED or GTF and do not select markers again.

`RunHspeBatch` uses one CPU and fits up to `hspe_batch_size` samples (default
100). The default worker memory is 4 GB. The existing `hspe_cpu`, `hspe_memory`,
and `hspe_disk_gb` settings apply to preparation and gene-type filtering. Increasing `hspe_cpu` does
not increase worker concurrency. Terra/Cromwell scheduling and cloud quotas
control how many batches run together. A failed batch can retry independently.

Each sample seed is derived from the workflow seed and its UTF-8 sample ID.
It does not depend on batch size, sample order, or worker scheduling. These seeds
differ from the previous single-call random sequence, so old results need not
be bitwise identical. The HSPE loss, optimizer, and stopping rules are unchanged.

`MergeHspeBatches` rejects duplicate, missing, or unexpected samples, mismatched
cell-type columns, invalid proportions, inconsistent seeds, and version conflicts.
It restores the original BED sample order. Its `hspe_sample_diagnostics` output
contains each sample ID, seed, iteration count, convergence code, and optimized
loss. Code 0 means convergence; code 1 means the iteration limit was reached.
The `hspe_log` output combines preparation, worker, and merge logs. Metadata
records batch size, batch count, selected-marker count, and seed strategy.
Batch RDS files are internal task outputs, not new public workflow outputs.

Major-lineage grouping and the cohort-mean abundance filter run after merging.
TCA still fits the full cohort on linear CPM. Supplying precomputed proportions
skips all three HSPE stages. The GitHub Actions fixture uses 12 samples in
batches of five, which also tests the final partial batch.

## Public input contracts

The tables below use the public WDL input names. `None` means that an optional
input has no default file. `Required` means that the workflow call must supply
a value.

### Standalone deconvolution inputs

| Input | Type | Default | Contract |
| --- | --- | --- | --- |
| `expression` | `File` | Required | Linear-CPM BED described below. |
| `gtf` | `File` | Required | Gene annotation used to map gene IDs to symbols in hspe mode. The WDL requires it in both modes. |
| `lm22` | `File` | Required | LM22 reference matrix. The workflow requires it in both proportion modes. |
| `precomputed_proportions` | `File?` | `None` | Optional sample-by-LM22 proportion matrix. If provided, the workflow skips hspe. |
| `covariates` | `File?` | `None` | Optional TCA covariates with `sample_id` first. Samples must match the expression order. Values must be finite numeric values, with no intercept or constant column. |
| `deconvolution_docker_image` | `String` | `"ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main"` | Container for all deconvolution tasks. |
| `preemptible_attempts` | `Int` | `2` | Global preemptible-attempt value for all tasks. |
| `max_retries` | `Int` | `2` | Global retry value for all tasks. |
| `min_lm22_overlap` | `Float` | `0.80` | Required LM22 gene-overlap fraction in `(0, 1]`. |
| `hspe_marker_fraction` | `Float` | `0.10` | hspe marker fraction in `(0, 1]`. |
| `hspe_quantile_normalize` | `Boolean` | `false` | If `true`, quantile-normalize the joined LM22 and bulk profiles before hspe. |
| `hspe_batch_size` | `Int` | `100` | Maximum samples per HSPE batch; must be positive. |
| `hspe_batch_memory` | `String` | `"4 GB"` | Memory for each single-CPU HSPE worker. |
| `hspe_batch_disk_gb` | `Int` | `10` | Local disk in GB for each HSPE worker. |
| `group_mean_threshold` | `Float` | `0.0001` | Finite nonnegative cohort-mean threshold for major-group retention. |
| `zero_floor` | `Float` | `0.000001` | Finite positive replacement for exact zero values in retained proportions. |
| `tca_max_iters` | `Int` | `10` | Positive TCA iteration limit. |
| `tca_parallel` | `Boolean` | `false` | If `true`, enable TCA parallel execution for model fitting and tensor export. CPU allocation is independent. |
| `random_seed` | `Int` | `20260901` | Positive random seed for HSPE and TCA. |
| `log2_pseudocount` | `Float` | `0.0` | Finite nonnegative value added to bulk CPM and LM22 before their hspe-only log2 transforms. |
| `gene_type` | `Array[String]` | `["protein_coding", "lncRNA"]` | Exact GTF gene types to retain before hspe and TCA. Values must be non-empty and unique. |
| `hspe_cpu` | `Int` | `4` | CPU count for the hspe task. |
| `hspe_memory` | `String` | `"32 GB"` | Memory for the hspe task. |
| `hspe_disk_gb` | `Int` | `100` | Local disk in GB for the hspe task. |
| `proportions_cpu` | `Int` | `2` | CPU count for proportion processing. |
| `proportions_memory` | `String` | `"16 GB"` | Memory for proportion processing. |
| `proportions_disk_gb` | `Int` | `50` | Local disk in GB for proportion processing. |
| `fit_cpu` | `Int` | `16` | CPU count for TCA fitting. |
| `fit_memory` | `String` | `"256 GB"` | Memory for TCA fitting. |
| `fit_disk_gb` | `Int` | `750` | Local disk in GB for TCA fitting. |
| `export_cpu` | `Int` | `8` | CPU count for TCA BED export. |
| `export_memory` | `String` | `"256 GB"` | Memory for TCA BED export. |
| `export_disk_gb` | `Int` | `500` | Local disk in GB for TCA BED export. |
| `gene_summary_memory` | `String` | `"8 GB"` | Memory for the gene summary task. It uses one CPU and `export_disk_gb` for disk space. |
| `manifest_cpu` | `Int` | `4` | CPU count for the deconvolution manifest. |
| `manifest_memory` | `String` | `"32 GB"` | Memory for the deconvolution manifest. |
| `manifest_disk_gb` | `Int` | `100` | Local disk in GB for the deconvolution manifest. |

### Integrated cell-type eQTL inputs

| Input | Type | Default | Contract |
| --- | --- | --- | --- |
| `expression` | `File` | Required | Linear-CPM BED described below. |
| `gtf` | `File` | Required | Gene annotation used to map gene IDs to symbols in hspe mode. The WDL requires it in both modes. |
| `lm22` | `File` | Required | LM22 reference matrix. The workflow requires it in both proportion modes. |
| `precomputed_proportions` | `File?` | `None` | Optional sample-by-LM22 proportion matrix. If provided, the workflow skips hspe. |
| `deconvolution_covariates` | `File?` | `None` | Optional TCA covariates. This is the integrated alias of standalone `covariates`. |
| `SampleList` | `File` | Required | Sample list passed to every scattered expression-QTL call. |
| `AdditionalCovariates` | `File` | Required | QTL covariates merged with selected phenotype PCs in each branch. |
| `OutputPrefix` | `String` | Required | Basename-safe token. It must start with an ASCII letter or number. It can contain only letters, numbers, dots, underscores, and hyphens. Each scatter call adds the cell-type slug. |
| `deconvolution_docker_image` | `String` | `"ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main"` | Container for deconvolution and integration tasks. |
| `qtl_docker_image` | `String` | `"ghcr.io/aou-multiomics-analysis/prepare_qtl:main"` | Container for scattered expression-QTL tasks. |
| `preemptible_attempts` | `Int` | `2` | Global preemptible-attempt value for deconvolution and QTL tasks. |
| `max_retries` | `Int` | `2` | Global retry value for deconvolution and QTL tasks. |
| `min_lm22_overlap` | `Float` | `0.80` | Required LM22 gene-overlap fraction in `(0, 1]`. |
| `hspe_marker_fraction` | `Float` | `0.10` | hspe marker fraction in `(0, 1]`. |
| `hspe_quantile_normalize` | `Boolean` | `false` | If `true`, quantile-normalize the joined LM22 and bulk profiles before hspe. |
| `hspe_batch_size` | `Int` | `100` | Maximum samples per HSPE batch; must be positive. |
| `hspe_batch_memory` | `String` | `"4 GB"` | Memory for each single-CPU HSPE worker. |
| `hspe_batch_disk_gb` | `Int` | `10` | Local disk in GB for each HSPE worker. |
| `group_mean_threshold` | `Float` | `0.0001` | Finite nonnegative cohort-mean threshold for major-group retention. |
| `zero_floor` | `Float` | `0.000001` | Finite positive replacement for exact zero values in retained proportions. |
| `tca_max_iters` | `Int` | `10` | Positive TCA iteration limit. |
| `tca_parallel` | `Boolean` | `false` | If `true`, enable TCA parallel execution for model fitting and tensor export. CPU allocation is independent. |
| `random_seed` | `Int` | `20260901` | Positive random seed for HSPE and TCA. |
| `log2_pseudocount` | `Float` | `0.0` | Finite nonnegative value added to bulk CPM and LM22 before their hspe-only log2 transforms. |
| `gene_type` | `Array[String]` | `["protein_coding", "lncRNA"]` | Exact GTF gene types to retain before hspe and TCA. Values must be non-empty and unique. |
| `hspe_cpu` | `Int` | `4` | CPU count for the hspe task. |
| `hspe_memory` | `String` | `"32 GB"` | Memory for the hspe task. |
| `hspe_disk_gb` | `Int` | `100` | Local disk in GB for the hspe task. |
| `proportions_cpu` | `Int` | `2` | CPU count for proportion processing. |
| `proportions_memory` | `String` | `"16 GB"` | Memory for proportion processing. |
| `proportions_disk_gb` | `Int` | `50` | Local disk in GB for proportion processing. |
| `fit_cpu` | `Int` | `16` | CPU count for TCA fitting. |
| `fit_memory` | `String` | `"256 GB"` | Memory for TCA fitting. |
| `fit_disk_gb` | `Int` | `750` | Local disk in GB for TCA fitting. |
| `export_cpu` | `Int` | `8` | CPU count for TCA BED export. |
| `export_memory` | `String` | `"256 GB"` | Memory for TCA BED export. |
| `export_disk_gb` | `Int` | `500` | Local disk in GB for TCA BED export. |
| `gene_summary_memory` | `String` | `"8 GB"` | Memory for the gene summary task. It uses one CPU and `export_disk_gb` for disk space. |
| `manifest_cpu` | `Int` | `4` | CPU count for the deconvolution manifest. |
| `manifest_memory` | `String` | `"32 GB"` | Memory for the deconvolution manifest. |
| `manifest_disk_gb` | `Int` | `100` | Local disk in GB for the deconvolution manifest. |
| `scatter_cpu` | `Int` | `1` | CPU count for scatter validation and the QTL manifest. |
| `scatter_memory` | `String` | `"4 GB"` | Memory for scatter validation and the QTL manifest. |
| `scatter_disk_gb` | `Int` | `20` | Local disk in GB for scatter validation and the QTL manifest. |
| `eqtl_cpu` | `Int` | `8` | CPU count for each scattered expression-QTL call. |
| `eqtl_memory` | `Int` | `64` | Memory in GB for each scattered expression-QTL call. |
| `eqtl_disk_gb` | `Int` | `200` | Local disk in GB for each scattered expression-QTL call. |

The standalone workflow uses `covariates` for the optional TCA model matrix.
The integrated workflow passes `deconvolution_covariates` to that input. It
keeps this matrix separate from required `AdditionalCovariates`. The integrated
workflow also has two image inputs because the deconvolution and QTL tasks use
different dependency sets.

Both workflows use one global `preemptible_attempts = 2` setting and one global
`max_retries = 2` setting. The CPU, memory, and disk inputs are task-specific
runtime settings.

## Expression input

The `expression` input is a BED file with gene IDs and sample values:

```text
#chr  start  end  gene_id  sample_1  sample_2
```

The first four columns must be `#chr`, `start`, `end`, and `gene_id` in that
order. The file must have at least one sample column. Chromosome values, gene
IDs, and sample IDs must be non-empty. Gene IDs and sample IDs must be unique.
`start` must be a nonnegative integer that is less than `end`. The sample
values must be linear CPM values. They must be finite and nonnegative.

The `gene_type` input defaults to `protein_coding` and `lncRNA`. The workflow
uses gene-type filtering before hspe and TCA. Values must match the GTF `gene_type`
or `gene_biotype` attributes exactly. The list must contain at least one unique,
non-empty value. A requested type can be absent from the GTF, but the workflow
fails if no expression genes remain. The `gene_type_filter_report` output shows
the annotation and retention status for each input gene. The
`filtered_expression` output keeps the original coordinates and order for all
retained genes.

The workflow uses the GTF to map retained gene IDs to gene symbols for the LM22
intersection. It removes constant genes before TCA. For each modeled,
nonconstant gene, the exported cell-type BED preserves the input coordinate and
modeled-gene order.

`log2_pseudocount` defaults to `0` and applies only to hspe. In hspe
mode, zero bulk CPM requires a positive pseudocount. Negative CPM values are
always invalid. The pseudocount must be finite and nonnegative. TCA accepts
zero CPM and does not use this pseudocount.

For proportion estimation, hspe uses `log2(CPM + log2_pseudocount)` for the
bulk matrix and `log2(LM22 + log2_pseudocount)` for the reference matrix.
Duplicate gene symbols are combined in linear CPM before the bulk transform.
TCA uses valid, nonconstant gene-ID rows directly in linear CPM. Each exported
cell-type BED is on the linear-CPM model scale. The workflow does not store a
full transformed bulk-expression matrix. Marker-only batches are temporary
HSPE task files. The export task rebuilds the
linear TCA view from `filtered_expression` and removes the same constant genes
before tensor extraction.

For compatibility, the integrated workflow sends each linear-scale TCA BED to
the `Log2CpmBed` input of `prepare_eQTL.wdl`. That input skips count
normalization and transformation. The QTL step applies INT to one branch and
direct centering and scaling to the other branch.

## Proportion modes

LM22 is required. If `precomputed_proportions` is not provided, the workflow
uses LM22 to estimate the 22 cell-type proportions with hspe. If
`precomputed_proportions` is provided, the workflow skips hspe and uses the
provided precomputed LM22 proportions.

- Provide only `lm22` to estimate proportions with hspe.
- Provide both inputs to use precomputed proportions and skip hspe.

The workflow accepts the official LM22 first-column header, `Gene symbol`. It
also accepts the canonical header, `gene_symbol`, and standardizes either form
to `gene_symbol` inside the hspe task.

### Precomputed LM22 schema

The precomputed file must have `sample_id` as the first column. Sample IDs must
be unique and non-empty. The values must be finite and nonnegative. Each row
must sum to one within `1e-6`. The other columns must be exactly the 22 standard
LM22 columns listed below. Their order in the input file can differ. The sample
row order must be exactly the same as the expression BED sample-column order.
The workflow does not reorder samples.

```text
B cells naive
B cells memory
Plasma cells
T cells CD8
T cells CD4 naive
T cells CD4 memory resting
T cells CD4 memory activated
T cells follicular helper
T cells regulatory (Tregs)
T cells gamma delta
NK cells resting
NK cells activated
Monocytes
Macrophages M0
Macrophages M1
Macrophages M2
Dendritic cells resting
Dendritic cells activated
Mast cells resting
Mast cells activated
Eosinophils
Neutrophils
```

Both modes combine LM22 values into major lineages. The major lineages include
separate CD4 T-cell and CD8 T-cell groups. The workflow removes a group when
its cohort mean is less than `group_mean_threshold`. No cell type is forced to
remain. If fewer than two major groups pass `group_mean_threshold`, the
workflow fails. `zero_floor` replaces exact zero values only in retained
proportions. It does not clamp other small positive values. The workflow then
normalizes the TCA weights by sample.

This example runs the standalone workflow with hspe:

```json
{
  "CellTypeDeconvolution.expression": "gs://bucket/expression.cpm.bed.gz",
  "CellTypeDeconvolution.gtf": "gs://bucket/genes.gtf.gz",
  "CellTypeDeconvolution.lm22": "gs://bucket/LM22.txt",
  "CellTypeDeconvolution.log2_pseudocount": 0.0
}
```

This example runs the standalone workflow with precomputed proportions:

```json
{
  "CellTypeDeconvolution.expression": "gs://bucket/expression.cpm.bed.gz",
  "CellTypeDeconvolution.gtf": "gs://bucket/genes.gtf.gz",
  "CellTypeDeconvolution.lm22": "gs://bucket/LM22.txt",
  "CellTypeDeconvolution.precomputed_proportions": "gs://bucket/lm22-proportions.tsv",
  "CellTypeDeconvolution.log2_pseudocount": 0.0
}
```

For the integrated workflow, use the `PrepareCellTypeEqtlWorkflow` input
prefix. Add the required `SampleList`, `AdditionalCovariates`, and
`OutputPrefix` inputs. For example:

```json
{
  "PrepareCellTypeEqtlWorkflow.expression": "gs://bucket/expression.cpm.bed.gz",
  "PrepareCellTypeEqtlWorkflow.gtf": "gs://bucket/genes.gtf.gz",
  "PrepareCellTypeEqtlWorkflow.lm22": "gs://bucket/LM22.txt",
  "PrepareCellTypeEqtlWorkflow.SampleList": "gs://bucket/samples.txt",
  "PrepareCellTypeEqtlWorkflow.AdditionalCovariates": "gs://bucket/covariates.tsv",
  "PrepareCellTypeEqtlWorkflow.OutputPrefix": "cohort",
  "PrepareCellTypeEqtlWorkflow.log2_pseudocount": 0.0
}
```

To use precomputed proportions in the integrated workflow, keep the required
`lm22` entry and add this entry:

```json
{
  "PrepareCellTypeEqtlWorkflow.precomputed_proportions": "gs://bucket/lm22-proportions.tsv"
}
```

`AdditionalCovariates` is required in the integrated workflow. The workflow
merges it with the selected phenotype PCs for each output branch.
`deconvolution_covariates` is optional and has a different role in the TCA
model.

## Per-gene CPM summaries

Both entry points publish `cell_type_gene_summary` as
`cell_type_gene_summary.tsv.gz`, and `gene_summary_log`. The
`SummarizeCellTypeBeds` task runs after TCA BED export in both HSPE and
precomputed-proportion modes. It matches localized BED basenames to the
inventory and reads one BED at a time in blocks of 256 genes. The task uses
one CPU, `gene_summary_memory` (default `"8 GB"`), `export_disk_gb`, and the
global retry settings. This adds one task; it does not add a scatter or
change TCA fitting.

The table has one row per cell type and gene. Cell types follow inventory
order; genes follow BED order. Gene IDs retain their version suffixes. The
columns are:

- `cell_type`, `n_samples`, `scale` (always `cpm`).
- `#chr`, `start`, `end`, `gene_id` from the BED.
- `mean_cpm`, `median_cpm` across samples.
- `sd_cpm`: sample standard deviation, with denominator `n_samples - 1`.
- `se_mean_cpm`: `sd_cpm / sqrt(n_samples)`.
- `q1_cpm`, `q3_cpm`: 25th and 75th percentiles, using R quantile type 7.
- `iqr_cpm`: `q3_cpm - q1_cpm`.

The task summarizes TCA-estimated linear CPM before downstream sample
selection, connectivity-outlier removal, INT, or scaling. It includes all
samples and genes present in the exported BEDs. It does not change the BEDs,
apply a log transform, remove zeros, or clamp negative TCA estimates. Invalid
or missing values cause an error instead of being silently removed. For one
sample, SD and standard error are `NA`; for a constant gene with two or more
samples, both are zero.

The standard error describes between-sample variation in the estimated
values, under an independent-sample assumption. It does not include TCA model
uncertainty and is not a standard error for the median. Related or repeated
samples need a separate dependence-aware analysis. These are CPM summaries,
not TPM summaries. Before comparing external datasets, check expression
units, gene identifiers, cell-type definitions, and processing differences.

The summary is a separate workflow output. The existing BED inventory and
QTL manifest keep their current schemas. GitHub Actions checks the summary
against the exported BEDs in both end-to-end smoke runs.

## Independent filtering and outputs

The workflow does independent connectivity filtering for each cell type and
for each INT and scaled branch. Thus, a cell type can have different retained
samples in its INT and scaled outputs. Two cell types can also have different
retained samples. Review the matching connectivity-outlier report before QTL
analysis.

The integrated workflow publishes ten ordered file arrays:

- `int_beds`
- `scaled_beds`
- `int_phenotype_pcs`
- `int_phenotype_pcs_all`
- `scaled_phenotype_pcs`
- `scaled_phenotype_pcs_all`
- `int_merged_covariates`
- `scaled_merged_covariates`
- `int_connectivity_outliers`
- `scaled_connectivity_outliers`

These `Array[File]` outputs are authoritative. Their order matches `cell_types`,
`cell_type_slugs`, and the manifest rows. The manifest stores stable basenames
because a WDL task cannot know the final Terra or Cromwell output URI. Use the
arrays to get the localized or cloud file paths.

The manifest does not contain raw or residualized BED files. The integrated
workflow exposes the deconvolution proportions, TCA model, weights, filter
report, BED inventory, reconstruction results, QC results, parameters, logs,
and provenance outputs separately.

The export task exposes its dynamically generated BED files with a WDL
`glob()`. The BED array order is not authoritative. The integration stage
matches each BED basename to the ordered inventory before it scatters the QTL
preparation tasks. The export stage also writes each BED SHA-256 checksum to
the inventory while the file is local. The manifest task uses these recorded
checksums and does not localize the large BED files onto a second worker.

## Terra file handling and validation

Optional HSPE metadata and TCA covariates remain `File?` values until command
rendering. Do not convert them to WDL `String` declarations: a string can
retain a `gs://` URI after Cromwell localizes the file. The command uses the
localized optional file directly, escapes apostrophes for the shell, and
omits the corresponding argument when the file is absent. The manifest task
checks that supplied HSPE metadata is readable before starting R.

Files must be created in task scope, not workflow scope. The static check
`python scripts/check_wdl_file_scope.py workflows` inspects all local WDL
sources. It rejects file-writing functions in workflow input defaults,
declarations, call inputs, outputs, scatters, and conditionals. Task-level
file creation is permitted. GitHub Actions runs this check and tests that
model file localization before command rendering, with optional files both
present and absent.

These local and CI checks do not execute Terra-managed Cromwell. The complete
workflow with this localization fix has not been tested on Terra. A Terra
submission requires user approval.

## Branch integration and repository ownership

Merge `feat/log2-cpm-bed-input` first. Merge
`codex/cell-type-specific-expression` after that branch reaches `main`. The
cell-type branch depends on the pre-normalized expression interface.

The old CellTypeDeconvolution repository is deprecated only after these
canonical workflows reach `prepare_QTL` on `main`. Do not delete the old
repository. Add its deprecation notice only after this integration is on
`main`.

## QTL manifest

`cell_type_qtl_manifest` is a tab-separated file. It has one row for each
retained cell type and these exact columns:

| Column | Meaning |
| --- | --- |
| `cell_type` | Cell-type display name. |
| `cell_type_slug` | Stable identifier used in filenames. |
| `int_bed` | Rank-based inverse-normal-transformed expression BED. |
| `scaled_bed` | Centered and scaled expression BED. |
| `int_phenotype_pcs` | Selected INT phenotype PCs. |
| `int_phenotype_pcs_all` | Complete INT phenotype PCs. |
| `scaled_phenotype_pcs` | Selected scaled phenotype PCs. |
| `scaled_phenotype_pcs_all` | Complete scaled phenotype PCs. |
| `int_merged_covariates` | Additional covariates merged with selected INT PCs. |
| `scaled_merged_covariates` | Additional covariates merged with selected scaled PCs. |
| `int_connectivity_outliers` | Connectivity-outlier report for the INT branch. |
| `scaled_connectivity_outliers` | Connectivity-outlier report for the scaled branch. |
