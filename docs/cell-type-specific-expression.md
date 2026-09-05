# Cell-type-specific expression workflows

These WDL 1.0 workflows prepare cell-type-specific expression from human
whole blood. `prepare_QTL` is the only maintained and canonical source for
this implementation.

## Start here: required files

Use this checklist for the integrated human whole-blood pipeline,
`PrepareCellTypeEqtlWorkflow`. Input names are case-sensitive.

| Input | Required? | File contents and scale | Use |
| --- | --- | --- | --- |
| `expression` | Yes | Gene-by-sample BED; finite, nonnegative **linear CPM**, not counts, TPM, log2 CPM, INT, or scaled values. | Supplies expression for HSPE and TCA. |
| `gtf` | Yes, in both proportion modes | GTF with gene records, `gene_id`, and `gene_type` or `gene_biotype`. HSPE also needs usable `gene_name` values. | Filters gene types and maps gene IDs to LM22 symbols. |
| `lm22` | Yes, even with supplied proportions | Original positive, linear LM22 reference values; genes in rows and the 22 LM22 cell types in columns. Do not pre-log this file. | Supplies the HSPE reference when proportions are estimated. |
| `SampleList` | Yes | One sample ID per line, with no header or an optional `research_id`, `sample_id`, or `ID` header. | Selects and orders samples for downstream eQTL preparation, **not** for HSPE or TCA. |
| `AdditionalCovariates` | Yes | Sample-by-covariate TSV with `sample_id`; numeric covariate columns. | Merges supplied covariates with phenotype PCs for each cell type and output branch. |
| `precomputed_proportions` | No | Sample-by-cell-type TSV; `sample_id` first, followed by all 22 LM22 columns. Fractions, not percentages. | Skips HSPE estimation. Grouping and TCA still run. |
| `deconvolution_covariates` | No | Sample-by-covariate TSV; `sample_id` first, followed by finite numeric columns. | Supplies covariates to the TCA model. It is not a substitute for `AdditionalCovariates`. |
| `haemopedia_counts` | No | Human Haemopedia **raw counts**, gene ID in the first column and sorted reference samples in the remaining columns; TSV or TSV.gz. | Adds reference expression filtering and comparison after TCA export. Do not supply CPM, TPM, or log values here. |

Also supply `OutputPrefix`, a safe filename prefix such as `whole_blood`.
The other integrated inputs have defaults. See the
[complete input tables](#public-input-contracts) for parameters and resources.

For **standalone deconvolution**, supply `expression`, `gtf`, and `lm22`.
Optional TCA covariates are named `covariates` in that workflow. It does not
take `SampleList`, `AdditionalCovariates`, or `OutputPrefix`.

For **eQTL preparation of one existing matrix**, use `eQTLPrepareData`.
Its inputs differ from the integrated pipeline. See
[CountGCT, CpmBed, and Log2CpmBed](#choose-the-expression-input-for-eqtl-preparation).

Quick links:

- [Scales and transformations](#scales-and-transformations)
- [File formats and sample matching](#file-formats-and-sample-matching)
- [Expression validation and gene filtering](#expression-input)
- [Proportions and major cell types](#proportion-modes)
- [Terra input examples](#terra-input-examples)
- [Output selection](#which-output-should-i-use)
- [Terra validation limits](#terra-file-handling-and-validation)

## Scales and transformations

Supply linear CPM to the deconvolution workflow. The workflow creates the
log2 values needed by HSPE internally. It does **not** pass those log2 values
to TCA.

```text
Gene-ID BED: linear CPM
  → GTF gene-type filter → filtered_expression: linear CPM
      ├─ HSPE: gene symbols → log2(CPM + p) → LM22 marker batches → proportions
      └─ TCA: linear CPM + grouped proportions → fit → model cleanup → cell-type BEDs: linear CPM
                                                  ├─ CPM gene summaries
                                                  └─ per-cell negative/reference filter → CpmBed → eQTL preparation
                                                       ├─ ranks → INT BED
                                                       └─ log2(CPM + 1) → center/scale → scaled BED
```

The HSPE branch is skipped when `precomputed_proportions` is supplied.
`p` means `log2_pseudocount`, whose default is `0`.

| Stage or file | Input scale | Handling | Output scale |
| --- | --- | --- | --- |
| Gene-type filtering | Linear CPM | Retains selected GTF gene types; preserves expression values. | Linear CPM |
| HSPE bulk input | Linear CPM | Sums gene IDs that map to the same symbol, then applies `log2(CPM + p)`. Uses shared LM22 genes and selected markers. | Log2 expression used internally |
| HSPE reference input | Original linear LM22 values | Applies `log2(LM22 + p)` internally. LM22 is a signature matrix, not this cohort's CPM matrix. | Log2 reference values used internally |
| HSPE estimates or supplied proportions | Fractions | Combines LM22 subtypes, filters major groups, floors exact zeros in retained groups, and normalizes each sample's weights to sum to one. | Fractions |
| TCA fitting and export | Linear CPM and grouped fractions | Removes constant input genes, fits the cohort model, and excludes numerically singular genes before extraction. Does not apply the HSPE log transform. | Estimated cell-type expression on the linear-CPM model scale |
| Gene summary | Exported TCA CPM | Computes mean, median, SD, standard error of the mean, quartiles, and IQR. | CPM units, except identifiers and sample counts |
| Reference comparison | Haemopedia raw counts and exported TCA CPM | TMM-normalizes the full reference matrix, calculates linear CPM, then compares sample means of `log2(CPM + 1)`. | Mean log2(CPM + 1) for comparisons; filtered BEDs retain linear CPM |
| eQTL INT branch | Supplied expression values | Applies rank-based inverse normal transformation per gene across selected samples. | Dimensionless INT values |
| eQTL scaled branch from `CpmBed` | Linear CPM | Applies `log2(CPM + 1)`, then centers and scales each gene across selected samples. | Dimensionless standardized values |

HSPE quantile normalization is **off by default**. If enabled, it applies to
the joined log2 bulk and reference profiles, before marker selection and
batching. It does not change the linear CPM supplied to TCA.

The two pseudocount uses are separate:

- `log2_pseudocount` affects only the HSPE bulk and reference transformations.
  It does not change TCA expression or exported TCA BEDs.
- The eQTL `CpmBed` scaled branch uses a fixed `+1`. It does not use
  `log2_pseudocount`. `log2(CPM + 1)` is not the same as `log2(CPM) + 1`.

At the default `p = 0`, HSPE checks the retained bulk CPM matrix for zeros
before the LM22 intersection. A zero outside the eventual marker set can
therefore stop HSPE preparation. Use a positive `p` when that matrix contains
zeros. LM22 values must be strictly positive even when `p` is positive.
TCA itself accepts zero CPM. Neither stage accepts negative bulk CPM input.

There is no automatic scale detection. Positive log2 values can look like
linear values to a numeric validator. Confirm the scale before submission;
do not choose an input name from the filename alone.

### Choose the expression input for eQTL preparation

The reusable `workflows/expression/prepare_eQTL.wdl` workflow accepts exactly
one of these three expression inputs:

| Input | Supplied values | Required annotation | Processing before scaled output |
| --- | --- | --- | --- |
| `CountGCT` | Raw gene counts in GCT format | `AnnotationGTF` required | Count filtering → TMM → CPM → `log2(CPM + 1)` → center and scale |
| `CpmBed` | Finite, nonnegative linear CPM in BED format | No GTF; supplying one is rejected | Skip count normalization; apply `log2(CPM + 1)` → center and scale |
| `Log2CpmBed` | Finite, already-log2 expression in BED format | No GTF; supplying one is rejected | Skip count normalization and log transformation; center and scale directly |

For counts, the current gene filter retains genes with counts greater than
6 in at least 20% of selected samples. BED modes do not repeat this filter.
Both BED modes preserve the supplied coordinates and reject zero-variance
genes. Finite negative log2 values are valid in `Log2CpmBed`; negative linear
CPM values are invalid in `CpmBed`.

Every mode creates an INT branch without a preceding log transform. The raw
output preserves supplied values for BED inputs; it contains TMM-normalized
CPM for count inputs. Thus, a file named `.raw.bed.gz` is not necessarily a
raw-count matrix. Its scale depends on the input mode.

The integrated pipeline automatically passes filtered TCA BEDs through `CpmBed`.
You do not set `CpmBed` or `Log2CpmBed` at its top level. Its public
`expression` input must still be linear CPM. In the reusable single-matrix
workflow, `SampleList`, `OutputPrefix`, `memory`, `disk_space`, and
`num_threads` are required; `AdditionalCovariates` is optional there.

## File formats and sample matching

Use actual tab characters for BED, TSV, and GTF fields. The examples below
show formats, not a complete biological test dataset. Localized files may
have different paths in Terra, but sample and gene identifiers must remain
unchanged.

### Expression BED

```text
#chr	start	end	gene_id	sample_1	sample_2	sample_3
chr1	999	1000	gene_A	1	7	3
chr2	1999	2000	gene_B	8	2	5
```

Rows are genes; sample columns contain linear CPM for `expression` or
`CpmBed`. Use the same layout, with already-transformed values, for
`Log2CpmBed`. BED coordinates use a zero-based start and an exclusive end.
The deconvolution workflow preserves these coordinates; it does not rebuild
them from the GTF. Use coordinates appropriate for the later QTL analysis.
Plain BED and gzip-compressed BED files are supported.

### Gene annotation GTF

```text
chr1	example	gene	1000	1500	.	+	.	gene_id "gene_A"; gene_name "ABCB4"; gene_type "protein_coding";
chr2	example	gene	2000	2400	.	+	.	gene_id "gene_B"; gene_name "EXAMPLE_LNC"; gene_type "lncRNA";
```

GTF records have nine tab-separated fields. The parser uses `gene` records,
not exon-only records. It accepts plain or gzip-compressed GTF files.
Use gene IDs that match the expression BED exactly, including version
suffixes. The pipeline does not remove version suffixes to make a match.
`gene_type` is preferred; `gene_biotype` is the fallback attribute.

The GTF is needed even with precomputed proportions because gene-type
filtering still runs. In HSPE mode, `gene_name` supplies the symbols used to
match LM22. Missing symbols prevent HSPE use of those genes, but do not by
themselves exclude an otherwise retained gene from TCA.

### LM22 reference

The first column is `Gene symbol` or `gene_symbol`. The remaining columns
must be all 22 standard LM22 cell types, with their original names, including
spaces and parentheses. Genes are rows. Values must be finite and strictly
positive; gene symbols must be unique. Column order can differ.

Supply the original linear reference values. Do not apply a log transform,
convert them to percentages, or treat them as the bulk CPM input. The
workflow performs its own log transformation. In precomputed-proportion
mode, `lm22` remains a required WDL input even though HSPE does not use it.

### Sample list

```text
sample_1
sample_2
sample_3
```

The header is optional. The first row can be `research_id`, `sample_id`, or
`ID`, in any letter case. A recognized header is removed only when it is not
an actual sample column in the expression file. The log records its removal.
Other first-row values are treated as sample IDs, not guessed headers.

Use exactly one column. Each ID must occur once and must match an expression
sample column. Blank IDs and empty lists cause an error. Sample order and
leading zeros are preserved. All three eQTL expression input modes require
every listed sample to be present; raw-count mode no longer silently drops
unmatched IDs.
The list selects and orders samples only after deconvolution. HSPE, TCA,
and the CPM gene summaries use all samples in the input expression BED.
To change that modeling cohort, prepare the expression BED and associated
proportion/covariate matrices with the intended samples before submission.

### QTL covariates: `AdditionalCovariates`

```text
sample_id	age	batch_B	genotype_pc1
sample_1	45	0	0.12
sample_2	52	1	-0.08
sample_3	39	0	0.03
```

Use one row per sample, a `sample_id` column, and uniquely named numeric
covariate columns. Encode categorical variables before submission; the
merge script does not encode them. Include every sample selected for QTL
preparation, with no missing covariate values. Avoid names such as `PC1`
that collide with the phenotype PC columns.

The merge joins `sample_id` to the phenotype-PC `ID` column. It is an inner
join: missing sample IDs can be dropped rather than cause an error.
Check coverage before submission. Input row order need not match the BED.
The merged output has covariates in rows, an `ID` first column, and sample
IDs in the remaining column headers. It sorts samples by ID. Always align
merged covariates and phenotype BEDs by sample ID, not by column position.

These covariates are not passed automatically to TCA. Merging them also does
not residualize the phenotype BEDs. The integrated workflow disables
residualized outputs; covariates are supplied separately for QTL analysis.

### Optional TCA covariates

Use the same sample-by-covariate layout shown above, with `sample_id` first.
The input is `deconvolution_covariates` in the integrated workflow and
`covariates` in standalone deconvolution. Unlike QTL covariates, its sample
rows must match **all expression BED sample columns in exactly the same
order**, not only the `SampleList` subset. Values must be finite numeric
values. Do not include an intercept or any constant covariate column.
The workflow does not log-transform or encode these covariates.

### Optional precomputed proportions

Use `sample_id` first, followed by all 22 LM22 column names listed in
[Precomputed LM22 schema](#precomputed-lm22-schema). Each row contains
fractions summing to one, within `1e-6`. For example, 15% is stored as
`0.15`, not `15`. Do not log-transform the values.

Include every expression sample in the same order as its BED column.
Do not supply only grouped CD4, CD8, or other major-lineage columns: the
workflow requires the 22 subtypes and performs the grouping itself.

### Raw-count GCT for the reusable eQTL workflow only

```text
#1.2
2	3
Name	Description	sample_1	sample_2	sample_3
gene_A	ABCB4	10	70	30
gene_B	EXAMPLE_LNC	80	20	50
```

The first two lines precede the table header; the second line gives the
gene and sample counts. The table must contain `Name`, `Description`, and
sample columns. Supply `AnnotationGTF` with matching gene IDs. This route
derives a one-base TSS interval from the GTF strand. Do not pass this GCT
to the deconvolution workflow's `expression` input.

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
| `gtf` | `File` | Required | Gene annotation for gene-type filtering in both modes and gene-ID-to-symbol mapping in hspe mode. |
| `lm22` | `File` | Required | LM22 reference matrix. The workflow requires it in both proportion modes. |
| `precomputed_proportions` | `File?` | `None` | Optional sample-by-LM22 proportion matrix. If provided, the workflow skips hspe. |
| `covariates` | `File?` | `None` | Optional TCA covariates with `sample_id` first. Samples must match the expression order. Values must be finite numeric values, with no intercept or constant column. |
| `haemopedia_counts` | `File?` | `None` | Raw human Haemopedia counts; absence means negative-only post-export filtering. |
| `reference_min_mean_log2_cpm1` | `Float` | `0.01` | Strict mean-log expression threshold in both datasets. |
| `reference_residual_cutoff` | `Float?` | `None` | Optional positive absolute standardized-residual cutoff; one pass, off by default. Requires the reference. |
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
| `gene_summary_memory` | `String` | `"8 GB"` | Memory for the gene summary, reference preparation, and BED filter tasks. Each uses one CPU and `export_disk_gb` for disk space. |
| `manifest_cpu` | `Int` | `4` | CPU count for the deconvolution manifest. |
| `manifest_memory` | `String` | `"32 GB"` | Memory for the deconvolution manifest. |
| `manifest_disk_gb` | `Int` | `100` | Local disk in GB for the deconvolution manifest. |

### Integrated cell-type eQTL inputs

| Input | Type | Default | Contract |
| --- | --- | --- | --- |
| `expression` | `File` | Required | Linear-CPM BED described below. |
| `gtf` | `File` | Required | Gene annotation for gene-type filtering in both modes and gene-ID-to-symbol mapping in hspe mode. |
| `lm22` | `File` | Required | LM22 reference matrix. The workflow requires it in both proportion modes. |
| `precomputed_proportions` | `File?` | `None` | Optional sample-by-LM22 proportion matrix. If provided, the workflow skips hspe. |
| `deconvolution_covariates` | `File?` | `None` | Optional TCA covariates. This is the integrated alias of standalone `covariates`. |
| `haemopedia_counts` | `File?` | `None` | Raw human Haemopedia counts; absence means negative-only post-export filtering. |
| `reference_min_mean_log2_cpm1` | `Float` | `0.01` | Strict mean-log expression threshold in both datasets. |
| `reference_residual_cutoff` | `Float?` | `None` | Optional positive absolute standardized-residual cutoff; one pass, off by default. Requires the reference. |
| `SampleList` | `File` | Required | Sample list with an optional `research_id`, `sample_id`, or `ID` header; passed to every scattered expression-QTL call. Does not subset HSPE or TCA. |
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
| `gene_summary_memory` | `String` | `"8 GB"` | Memory for the gene summary, reference preparation, and BED filter tasks. Each uses one CPU and `export_disk_gb` for disk space. |
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
nonconstant gene retained after post-fit numerical cleanup, the exported
cell-type BED preserves the input coordinate and modeled-gene order.

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
linear TCA view from `filtered_expression`, removes the same constant genes,
and applies the cleaned model's recorded numerical exclusions before tensor
extraction. See [post-fit cleanup](#post-fit-numerical-cleanup).

The integrated workflow sends each linear-scale TCA BED to the `CpmBed` input
of `prepare_eQTL.wdl`. This input skips count filtering, TMM normalization,
CPM calculation, and GTF mapping. The QTL step rank-normalizes the supplied
values for INT. For the scaled output, it applies `log2(CPM + 1)` before
centering and scaling. Connectivity-outlier checks, phenotype PCs, and merged
covariates remain separate for the two branches. TCA fitting, exported TCA
BED files, and gene summaries remain in linear CPM space.

The `Log2CpmBed` input remains available for files that are already in log2
space and never applies another log transform. Supply exactly one of
`CountGCT`, `CpmBed`, or `Log2CpmBed` to the reusable eQTL workflow.
`CpmBed` rejects negative or non-finite estimates before writing expression
outputs. It does not clamp negative TCA estimates to zero. If an exported TCA
BED contains negative values, eQTL preparation stops and reports example
gene IDs; the original TCA output remains unchanged.

Deploy the updated standard QTL image together with the updated WDL: older
images do not implement the `CpmBed` input. This change does not require a
different image for TCA. GitHub Actions checks the numerical transformation
from exported TCA BEDs to scaled QTL BEDs. The complete workflow with this
CPM-input change has not been tested on Terra; no Terra job was submitted.

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

The default `group_mean_threshold = 0.0001` means **0.01%**, not 1%.
It applies to each major group's mean fraction across the full expression
cohort, not to individual samples. Genes missing from LM22 are not removed
from TCA for that reason. HSPE alone uses the reference intersection; the
default `min_lm22_overlap = 0.80` is the fraction of genes in the supplied
LM22 reference that match the retained, symbol-mapped bulk matrix.

| Major group | LM22 columns that are summed |
| --- | --- |
| B cells | B cells naive; B cells memory; Plasma cells |
| CD4 T cells | T cells CD4 naive; T cells CD4 memory resting; T cells CD4 memory activated; T cells follicular helper; T cells regulatory (Tregs) |
| CD8 T cells | T cells CD8 |
| Gamma-delta T cells | T cells gamma delta |
| NK cells | NK cells resting; NK cells activated |
| Monocyte/myeloid | Monocytes; Macrophages M0; Macrophages M1; Macrophages M2 |
| Neutrophils | Neutrophils |
| Eosinophils | Eosinophils |
| Dendritic cells | Dendritic cells resting; Dendritic cells activated |
| Mast cells | Mast cells resting; Mast cells activated |

These are the pipeline's group definitions. The labels, particularly
Monocyte/myeloid, describe combined signatures rather than a measured count
of a single sorted population.

## Terra input examples

Replace `gs://bucket/...` with accessible files in your workspace storage.
The examples use `log2_pseudocount = 0.0`, so the retained bulk matrix must
have no zeros when HSPE runs. Set this parameter to a positive value if
needed. Do not add that pseudocount to the BED beforehand.

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

To run only eQTL preparation on an existing **linear CPM BED**, use this
input object with `workflows/expression/prepare_eQTL.wdl`:

```json
{
  "eQTLPrepareData.CpmBed": "gs://bucket/cd4_t_cells.bed.gz",
  "eQTLPrepareData.SampleList": "gs://bucket/samples.txt",
  "eQTLPrepareData.AdditionalCovariates": "gs://bucket/covariates.tsv",
  "eQTLPrepareData.OutputPrefix": "cohort.cd4_t_cells",
  "eQTLPrepareData.memory": 64,
  "eQTLPrepareData.disk_space": 200,
  "eQTLPrepareData.num_threads": 8
}
```

For an already-log2 BED, replace the `CpmBed` entry with `Log2CpmBed`.
Do not provide both. For counts, provide `CountGCT` and `AnnotationGTF`
instead of either BED input. `AdditionalCovariates` is optional only in this
single-matrix workflow. The integrated workflow requires it.

Use the WDL and images from a matching release. The integrated workflow's
`qtl_docker_image` must contain the `CpmBed` implementation added in
[PR #34](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/pull/34).
It uses `deconvolution_docker_image` for HSPE/TCA and `qtl_docker_image` for
eQTL preparation. A mutable `:main` tag can change; record the resolved
image digest and workflow revision for each run.

Memory input types differ. For example, `fit_memory` and `export_memory`
are strings such as `"256 GB"`. `eqtl_memory` is an integer such as `64`;
the expression task adds the GB unit. `tca_parallel = false` remains the
default even when multiple CPUs are allocated. Global `max_retries` does
not by itself request more memory on a retry.

Before a Terra submission, check:

1. The expression scale matches the selected workflow input.
2. BED and GTF gene IDs use matching identifiers and version suffixes.
3. Proportions and TCA covariates match all BED samples in the correct order.
4. QTL covariates cover every sample selected in `SampleList`.
5. The HSPE pseudocount is suitable for zeros in the retained bulk matrix.
6. The workflow revision and both image versions are recorded.

## Which output should I use?

| Output | Contents or scale | Use |
| --- | --- | --- |
| `tca_model` | Final cleaned TCA model RDS; gene-specific parameters for retained genes only. | Use for downstream tensor extraction or other analyses of the final model. |
| `tca_model_unfiltered` | Original fitted TCA model RDS, before numerical cleanup. | Audit the fit or repeat cleanup without fitting again. Do not substitute it for the final model. |
| `tca_numerical_excluded_genes` | Gene IDs, reasons, variance ranges, reciprocal condition numbers, and threshold. | Audit numerical exclusions after fitting. This is separate from the pre-fit constant-gene report, `tca_excluded_genes`. |
| `cell_type_beds` | One gene-by-sample BED per retained cell type; estimated linear CPM. | Compare cell-type expression profiles. Keep the matching gene IDs, sample IDs, and inventory. |
| `filtered_cell_type_beds` | Separate nonnegative, per-cell-type filtered BEDs; unchanged linear CPM values for retained rows. | These are the inputs to the integrated eQTL scatter. |
| `negative_expression_summary`, `reference_gene_comparison`, `reference_filter_metrics` | Negative-value counts, gene-level comparisons/exclusions, and per-cell filter/comparison statistics. | Audit the retained gene universe. See the reference filtering section below. |
| `cell_type_gene_summary` | Per-cell-type, per-gene CPM summaries across all exported samples. | Compare mean or median expression profiles; inspect variability. See the standard-error limitations below. |
| `estimated_proportions` | HSPE's sample-by-22-cell-type fractions; absent when precomputed proportions are used. | Inspect the proportion estimates before grouping. |
| `proportions_lm22`, `proportions_combined`, `tca_weights` | Original 22-type fractions, summed major groups, and filtered/renormalized model weights. | Inspect the changes made before TCA. These are not expression matrices. |
| `int_beds` | Dimensionless inverse-normal-transformed phenotypes after outlier removal. | Run QTL analysis with matching INT PCs, covariates, and sample IDs. |
| `scaled_beds` | Dimensionless standardized log2(CPM + 1) phenotypes after outlier removal. | Run QTL analysis with matching scaled PCs, covariates, and sample IDs. |
| `int_phenotype_pcs`, `scaled_phenotype_pcs` | Selected phenotype PC scores for the corresponding branch. | Inspect PCs and use the matching merged covariates. The `_all` outputs contain the full PC sets. |
| `int_merged_covariates`, `scaled_merged_covariates` | Supplied QTL covariates plus selected phenotype PCs; covariates in rows. | Adjust QTL association models. Align samples by ID. |
| `int_connectivity_outliers`, `scaled_connectivity_outliers` | Removed `SampleID` values and connectivity `Z_score`. | Check why sample sets differ across cell types or branches. |
| `cell_type_qtl_manifest` | One row per cell type, with a Terra entity ID, ten QTL output categories, source/filtered CPM BEDs and filter reports. | Import the TSV into a Terra data table or use the paths directly. |

Do not interpret INT values, scaled values, PCs, or proportions as CPM or
TPM. For external expression comparisons, start with `cell_type_beds` or
`cell_type_gene_summary`, not the transformed QTL BEDs. These outputs are
model estimates, not directly measured sorted-cell expression. Match gene
identifiers, cell-group definitions, and expression units before comparing
datasets. This pipeline does not calculate TPM or convert CPM to TPM.

The integrated workflow does not publish the nested eQTL `.raw` BEDs as
top-level outputs. Use its original TCA `cell_type_beds` for the untransformed
cell-type profiles. A separate run of `eQTLPrepareData` exposes `RawBedFile`;
that file has already been subset to its `SampleList`.

For troubleshooting, retain `gene_type_filter_report`, `cell_group_filter_report`,
`tca_excluded_genes`, HSPE overlap/marker diagnostics, reconstruction results,
QC plots, and task logs. TCA can export negative estimates; the original BEDs
and original summary keep them. The post-export filter removes a gene from the
affected cell type before its BED reaches `CpmBed`. Other eQTL input checks still
apply, so successful filtering does not guarantee every downstream call passes.

## Post-fit numerical cleanup

Both public workflows run `CleanTcaModel` after `FitTca` and before
`ExportTcaBeds`. The public `tca_model` output is
`tca_model_cleaned.rds`. Export and manifest generation use this final model.
The original `tca_model.rds` remains available as `tca_model_unfiltered`.
The cleanup task uses one CPU, 4 GB memory, and a 20 GB disk, with the global
preemptible and retry settings.

For each gene, cleanup squares the fitted standard deviations in
`sigmas_hat`. It then calculates:

```text
reciprocal_condition = min(cell-type variances) / max(cell-type variances)
exclude when reciprocal_condition < .Machine$double.eps
```

The threshold is approximately `2.220446049250313e-16`, the default
tolerance used by R's matrix solver. For the diagonal variance matrix used
by TCA, the ratio is its reciprocal condition number. An all-zero variance
row receives a ratio of zero and is excluded. Nonfinite or negative standard
deviations, overflow during squaring, and invalid model identifiers cause an
error instead of a silent exclusion. Cleanup also fails if no genes remain.
The shared residual standard deviation, `tau_hat`, must be finite and
positive because the unchanged TCA exporter divides by its square. An
invalid `tau_hat` stops cleanup; it does not trigger gene exclusions.

Cleanup does not refit TCA, change fitted values, add a pseudocount, alter
cell-type weights, or change the statistical model. It subsets all known
gene-specific parameter matrices, including covariate effects and p-values.
Sample-level matrices and the shared residual standard deviation remain
unchanged. The final RDS records the original gene order, excluded gene IDs,
filter method, and threshold in `gene_filter`.

Before extraction, the exporter checks the original variable-gene order
against that record. It removes only the recorded exclusions and aligns
the BED coordinates to the retained model rows. It does not use an unchecked
gene intersection. Cell-type BEDs, reconstruction QC, gene summaries, and
downstream QTL preparation therefore use the same retained gene set before
any QTL-specific processing. The earlier `filtered_expression` output is
still the GTF-filtered input BED, not a numerically cleaned BED.

`tca_numerical_excluded_genes.tsv` contains these columns:
`gene_id`, `reason`, `min_variance`, `max_variance`,
`reciprocal_condition`, and `threshold`. Variances are in squared CPM units;
the ratio and threshold are dimensionless. A run with no exclusions writes
a header-only report. `tca_cleanup_log` records the input, retained, and
excluded gene counts. Numerical filtering is not proof that the remaining
cell-specific estimates are biologically reliable, and it does not remove
negative expression estimates.

To clean an existing fit without refitting, run this command inside the
matching updated cell-type image:

```bash
Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/clean_tca_model.R \
  --model /path/to/tca_model.rds \
  --output-dir /path/to/cleaned
```

Use the resulting `tca_model_cleaned.rds` with the original GTF-filtered CPM
BED and matching sample weights when running the updated exporter. Do not
manually remove expression rows first: the exporter uses the recorded
original gene order to check alignment. Existing complete models without
cleanup metadata remain supported by the direct export script only when
their gene order exactly matches the variable expression rows. The public
workflows always run cleanup.

The cleanup change has local model and export regression tests. The complete
updated workflow has not been tested on Terra; no Terra job was submitted
for this change.

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

The summary is a separate workflow output. The original BED inventory keeps
its schema; the QTL manifest adds the filtering paths described below.
GitHub Actions checks the original summary
against the exported BEDs in both end-to-end smoke runs.

## Per-cell-type reference filtering

After export, every cell type is filtered independently. If any sample has a
CPM estimate below zero, remove that entire gene row from that cell type's
filtered BED. There is no tolerance for tiny negative values and no clipping.
The same gene can remain in other cell types. The negative summary records
counts, percentages, minimum CPM, and mean negative CPM before any removal.
Sample membership is unchanged; the check uses all exported samples, not only
the later eQTL `SampleList`.

Without `haemopedia_counts`, this is the only new filter. With the reference,
the workflow follows the supplied edgeR procedure on the full count matrix:

```r
y <- edgeR::DGEList(counts)
y <- edgeR::calcNormFactors(y)
reference_cpm <- edgeR::cpm(y, log = FALSE)
```

Only then are genes matched with the deconvolution output. Do not pre-normalize
the reference or restrict its gene universe before passing it as raw counts.
Join by gene ID; numeric version suffixes on Ensembl IDs can be removed for
matching, but original BED identifiers stay unchanged. Ambiguous duplicate
matching keys are errors. Missing reference genes are not measured zeros.

The workflow removes terminal replicate suffixes such as `.1`, then maps:

| Pipeline cell type | Reference populations |
| --- | --- |
| B cells | `NveB`, `MemB` |
| CD4 T cells | `CD4T` |
| CD8 T cells | `CD8T` |
| NK cells | `NK` |
| Monocyte/myeloid | `Mono`, `MonoNonClassical` |
| Neutrophils | `Neut` |
| Eosinophils | `Eo` |
| Dendritic cells | `myDC`, `myDC123`, `pDC` |

Each reference sample has equal weight, not each subtype. For example, five
naive-B samples contribute more than three memory-B samples. All expected
subpopulations must be present for a compared lineage. Mast and gamma-delta
T cells have no matching reference here: only their negative filter applies.
When a reference is supplied, it must contain every mapped lineage used by the
BED inventory. A missing required lineage is an input error; the workflow does
not silently switch that lineage to negative-only filtering. A reference can
contain only a subset of mapped lineages when the BED inventory uses only that
same subset. This requirement does not apply when the complete reference input
is omitted, or to mast and gamma-delta T cells.

For each eligible gene, compute **the mean of sample-level `log2(CPM + 1)`**
in each dataset. This differs from `log2(mean(CPM) + 1)` and from the existing
CPM-scale summary. Retain genes with mean log expression strictly greater than
`reference_min_mean_log2_cpm1` (default `0.01`) in both datasets. This is a
permissive expression threshold. Unmatched genes are excluded only in mapped
reference cell types. The filtered BEDs themselves remain linear CPM.

Compare deconvolution mean log expression (y) against reference mean log
expression (x) across genes. Reports include Pearson r, Spearman rho, OLS
R-squared, slope, intercept, gene counts, residuals and standardized residuals.
Plots include negative prevalence, a heatmap of up to 100 genes with the largest
negative fractions, and per-cell scatter and residual plots. Full gene results
are retained in the tables even when the heatmap shows only a subset.

`reference_residual_cutoff` is optional and **off by default**. If supplied,
remove genes whose absolute internally standardized OLS residual exceeds it.
The residual is `e / (s * sqrt(1 - h))`, where `s` is the residual standard error
and `h` is leverage. This is one pass using the baseline fit, not iterative
outlier removal. Baseline metrics stay separate from metrics on retained genes.
No cutoff is recommended as a universal biological threshold.

An insufficient or constant comparison gives an explicit unavailable status.
Requested residual removal fails when the residuals cannot be calculated.
A cell type with no remaining genes fails before eQTL preparation.

The B-cell model group includes plasma cells, absent from this reference. The
myeloid model group includes macrophage components, unlike the reference
monocyte populations. The dendritic group combines distinct subtypes. A large
residual can reflect these differences or real biology, not just estimation
error. Correlation across genes is descriptive, not donor-level validation.
R-squared here equals squared Pearson r and does not measure identity-line
agreement. Better correlation after residual filtering is selection-dependent.

Example additional Terra inputs:

```json
{
  "PrepareCellTypeEqtlWorkflow.haemopedia_counts": "gs://YOUR_BUCKET/GSE115736_Haemopedia-Human-RNASeq_raw.txt.gz",
  "PrepareCellTypeEqtlWorkflow.reference_min_mean_log2_cpm1": 0.01
}
```

The workflow publishes `filtered_cell_type_bed_inventory`,
`negative_expression_summary`, `reference_gene_comparison`,
`reference_filter_metrics`, `reference_filter_plots`, and `reference_filter_log`.
When a reference is supplied, `haemopedia_reference_summary`,
`haemopedia_reference_samples`, and `haemopedia_reference_metadata` document its
processing. Original BEDs, the TCA model, raw summaries and reconstruction QC
are preserved. The reference filter does not change TCA or HSPE calculations.

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
`cell_type_slugs`, and the manifest rows. The final QTL manifest preserves the
completed upstream calls' full output paths: `gs://` URLs on Terra and absolute
host paths in local smoke tests. The manifest task receives these paths as
`Array[String]` metadata, so Cromwell does not localize the files again. The
task validates the path format and array alignment; it does not read the file
contents. No files are copied or moved. No destination path input is needed.

The manifest includes original and filtered CPM BEDs, but not residualized BED files. The integrated
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

The maintained workflows are now in `prepare_QTL` on `main`. New users do
not need to merge the old development branches.

The expression-scale changes have this history:

- [Commit 9edfbdb, June 25, 2026](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/commit/9edfbdb43a8ad7c63f216b4224df0153639d8ad0)
  added `log2(CPM + 1)` before scaling in the count-input route. Immediately
  before that change, the scaled branch used linear CPM, not `log2(CPM)`.
- [Commit f8552c6, September 2, 2026](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/commit/f8552c6e7709ea02d4ae0acc02d5729e8af0084c)
  added `Log2CpmBed` and skipped the log transform for that input.
- [PR #34](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/pull/34)
  added `CpmBed` and routed linear TCA output through it. This corrected the
  earlier connection to `Log2CpmBed`, which had scaled TCA CPM directly.
  It did not change TCA fitting or export.

Historical integration order: `feat/log2-cpm-bed-input` first, then
`codex/cell-type-specific-expression` after the expression interface was
available. The old CellTypeDeconvolution repository was to be deprecated
only after the canonical workflows reached `prepare_QTL` on `main`.
Do not delete the old repository. Use `prepare_QTL` for current development;
the historical plan files are not current operator instructions.

## QTL manifest

`cell_type_qtl_manifest` is a tab-separated file. It has one row for each
retained cell type and these exact columns:

| Column | Meaning |
| --- | --- |
| `entity:cell_type_id` | First column; Terra row ID, equal to the cell-type slug (for example, `cd4_t_cells`). |
| `cell_type` | Cell-type display name. |
| `cell_type_slug` | Stable identifier used in filenames. |
| `int_bed` | Rank-based inverse-normal-transformed expression BED. |
| `scaled_bed` | Expression BED after log2(CPM + 1), centering, scaling, and connectivity-outlier removal. |
| `int_phenotype_pcs` | Selected INT phenotype PCs. |
| `int_phenotype_pcs_all` | Complete INT phenotype PCs. |
| `scaled_phenotype_pcs` | Selected scaled phenotype PCs. |
| `scaled_phenotype_pcs_all` | Complete scaled phenotype PCs. |
| `int_merged_covariates` | Additional covariates merged with selected INT PCs. |
| `scaled_merged_covariates` | Additional covariates merged with selected scaled PCs. |
| `int_connectivity_outliers` | Connectivity-outlier report for the INT branch. |
| `scaled_connectivity_outliers` | Connectivity-outlier report for the scaled branch. |
| `source_cpm_bed` | Original TCA BED for this cell type, before post-export filtering. |
| `filtered_cpm_bed` | Linear CPM BED supplied to eQTL preparation after per-cell filtering. |
| `negative_expression_summary` | Shared per-gene, per-cell-type negative-value report. |
| `reference_gene_comparison` | Shared comparison and gene-exclusion report. |
| `reference_filter_metrics` | Shared filter counts and comparison metrics. |

All file columns contain full paths, not basenames. Shared report paths repeat
across cell-type rows. Import this TSV with
Terra's **Import Data** action to create or update the `cell_type` table.
The pipeline does not import it automatically. See the
[Terra table format guide](https://support.terra.bio/hc/en-us/articles/6197368140955-How-to-make-a-data-table-from-scratch-or-a-template).

IDs are stable across runs. Importing another run with the same IDs into the
same table updates those rows. Use a separate table or change the IDs before
import if you need to keep multiple cohorts or runs in one workspace.
Keep the referenced submission files: the TSV contains links, not copies.
Call-cached outputs can link to an earlier run. The separate deconvolution
`output_manifest` is unchanged. No provenance file is added.
# Restart from a TCA model

Both entry points accept the optional `File? precomputed_tca_model` input.
For the complete workflow, set
`PrepareCellTypeEqtlWorkflow.precomputed_tca_model` to the cloud path of a model RDS.
It has priority over `precomputed_proportions`. If absent, the existing run mode
is unchanged.

Supply a pipeline model with `expression_scale = "cpm"`: either `tca_model`
(cleaned) or `tca_model_unfiltered` (fitted). Models on a log scale, or without
a recorded scale, are rejected. Numerical cleanup runs before export; an already
cleaned model keeps its recorded exclusions.

The expression BED must contain matching **linear CPM values**, the same samples
in the same order, and all retained model genes. Extra genes are ignored, and
gene rows and BED coordinates are aligned to the model. The workflow cannot
confirm that values came from the original fit; use the same dataset and scale.

Restart skips GTF gene filtering, HSPE, proportion processing, and TCA fitting.
The model's genes, cell types, weights, and covariates are authoritative. New
gene-type filters, proportion thresholds, fitting settings, and deconvolution
covariates do not change the saved model. The existing required GTF and LM22
inputs remain required for interface compatibility, but restart does not read
them for deconvolution. QTL sample and covariate inputs still apply downstream.

BED export, reference filtering, summaries, QTL preparation, and manifests still
run. Outputs from skipped stages are null. `tca_model` and `tca_weights` remain
available. QC marks unavailable fitting/proportion metrics as
`unavailable_model_restart`; the effective parameters do not claim new settings
were used to fit the supplied model.

GitHub Actions tests restart from the model produced by the existing smoke run.
It checks skipped stages, identical BED contents, and downstream QTL outputs.
This is not a complete Terra validation.
