# Cell-type-specific expression workflows

These WDL 1.0 workflows prepare cell-type-specific expression from human
whole blood. `prepare_QTL` is the only maintained and canonical source for
this implementation.

## Public workflows

Use
[`deconvolution.wdl`](../workflows/cell_type_specific_expression/deconvolution.wdl)
to run dtangle, process the proportions, fit TCA, and export one BED for each
retained major cell type.

Use
[`prepare_cell_type_eQTL.wdl`](../workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl)
to run the same deconvolution and then prepare QTL inputs for each retained
cell type. This workflow calls the existing single-matrix
[`prepare_eQTL.wdl`](../workflows/expression/prepare_eQTL.wdl). It does not
copy the expression-QTL tasks.

## Expression input

The `expression` input is a BED file with gene IDs and sample values:

```text
#chr  start  end  gene_id  sample_1  sample_2
```

The sample values must be linear CPM values. They must be finite and
nonnegative. The workflow does not filter the input to protein-coding genes.
It uses the GTF to map gene IDs to gene symbols for the LM22 intersection.

`log2_pseudocount` defaults to `0`. When `log2_pseudocount` is `0`, all CPM
values must be strictly positive. Zero CPM values are valid only when
`log2_pseudocount` is greater than zero. Negative CPM values are always
invalid. The pseudocount must be finite and nonnegative.

dtangle and TCA apply `log2(CPM + log2_pseudocount)` exactly once. Duplicate
gene symbols are first combined in linear CPM for the dtangle view. TCA uses
valid, nonconstant gene-ID rows. Each exported cell-type BED is in log2-CPM
space.

The integrated workflow sends each TCA BED to the `Log2CpmBed` input of
`prepare_eQTL.wdl`. The QTL step does not apply a second log2 transform. It
applies INT to one branch and direct centering and scaling to the other branch.

## Proportion modes

You must supply exactly one of lm22 or precomputed_proportions. Do not supply
both inputs.

- Supply `lm22` to estimate the 22 LM22 cell-type proportions with dtangle.
- Supply `precomputed_proportions` to use precomputed LM22 proportions and
  skip dtangle.

Both modes combine LM22 values into major lineages. The major lineages include
separate CD4 T-cell and CD8 T-cell groups. The workflow removes a group when
its cohort mean is less than `group_mean_threshold`. No cell type is forced to
remain. It applies `zero_floor` to retained proportions and then normalizes the
TCA weights by sample.

This example runs the standalone workflow with dtangle:

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

To use precomputed proportions in the integrated workflow, remove the `lm22`
entry and add this entry:

```json
{
  "PrepareCellTypeEqtlWorkflow.precomputed_proportions": "gs://bucket/lm22-proportions.tsv"
}
```

`AdditionalCovariates` is required in the integrated workflow. The workflow
merges it with the selected phenotype PCs for each output branch.
`deconvolution_covariates` is optional and has a different role in the TCA
model.

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
