# Cell-type-specific eQTL integration design

**Date:** 2026-09-02

**Status:** Approved in conversation; pending written-spec review

## Goal

Make `prepare_QTL` the only maintained source for the whole-blood dtangle and
TCA deconvolution pipeline. Add one end-to-end WDL 1.0 workflow that creates
cell-type-specific log2-CPM BED files, scatters the existing expression-QTL
preparation workflow over those files, and writes one portable manifest row per
retained cell type.

## Scope

This change will:

- Migrate the tested deconvolution implementation from
  `AoU-Multiomics-Analysis/CellTypeDeconvolution` commit `9a98d2e` into
  `prepare_QTL`.
- Preserve the validated dtangle, LM22 grouping, TCA, BED export, and QC
  behavior.
- Add a standalone deconvolution workflow in `prepare_QTL`.
- Add an end-to-end workflow that scatters expression-QTL preparation over the
  retained cell types.
- Add a manifest that links each cell type to its QTL-ready INT and scaled
  files.
- Add one user-configurable log2 pseudocount that defaults to zero.
- Add pinned container and GitHub Actions support for the migrated
  deconvolution code.
- Add a final deprecation notice to the old repository after the integrated
  workflow is green.

This change will not:

- Delete the old repository.
- Maintain synchronized copies of the deconvolution implementation.
- Add residualized BED files to the cell-type manifest.
- Change the statistical dtangle or TCA model.
- Reapply a log transform in expression-QTL preparation.
- Force any cell type to pass the existing cohort-mean retention threshold.

## Canonical repository layout

The migrated source will use these boundaries:

```text
workflows/cell_type_specific_expression/
  deconvolution.wdl
  prepare_cell_type_eQTL.wdl
  tasks/

scripts/cell_type_specific_expression/
  R/
  command-line scripts

envs/CellTypeSpecificExpression/
  Dockerfile
  environment.yml
```

Tests and fixtures will use matching cell-type-specific directories under
`tests/`. User documentation will be added to the repository workflow catalog
and a focused guide in `docs/`.

The existing `workflows/expression/prepare_eQTL.wdl` remains the canonical
single-matrix expression-QTL preparation workflow. The new end-to-end workflow
will import and call it. It will not copy its tasks.

## Public workflow entry points

### Standalone deconvolution

`workflows/cell_type_specific_expression/deconvolution.wdl` will expose the
current direct-CPM workflow after path and container-default updates. Its input
expression BED contains:

```text
#chr  start  end  gene_id  <sample columns...>
```

The sample values must be finite, nonnegative, linear CPM values. When
`log2_pseudocount` is `0`, all values must be strictly positive. Zero values are
valid only when `log2_pseudocount` is greater than zero. Negative values are
always invalid.

The standalone workflow will preserve both proportion modes:

- Estimate LM22 proportions with dtangle from the expression BED and GTF.
- Accept precomputed LM22 proportions instead of running dtangle.

It will preserve the current outputs, including proportions, retained-group
report, TCA model, one log2-CPM BED per retained group, BED inventory,
reconstruction QC, plots, logs, and provenance manifest.

### End-to-end cell-type eQTL preparation

`workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl` will call
the standalone deconvolution workflow, validate and label its BED outputs, and
scatter the existing `eQTLPrepareData` workflow over the retained groups.

The end-to-end workflow will require:

- The deconvolution expression BED.
- The GTF.
- Exactly one valid proportion mode: LM22 or precomputed proportions.
- A sample list for expression-QTL preparation.
- `AdditionalCovariates` for merging with phenotype PCs.
- An output prefix.

The workflow will expose:

```wdl
Float log2_pseudocount = 0.0
```

The value must be finite and greater than or equal to zero. The same value will
be used for the dtangle and TCA expression views and recorded in the effective
parameter output.

Optional deconvolution covariates remain separate from the required additional
QTL covariates. They have different statistical roles and will have distinct
input names.

## Data flow

1. The deconvolution subworkflow reads the linear-CPM BED.
2. dtangle maps gene IDs through the GTF, aggregates duplicate symbols in
   linear CPM, intersects the result with LM22, and applies
   `log2(CPM + log2_pseudocount)` exactly once.
3. TCA uses every valid, nonconstant gene-ID row and applies
   `log2(CPM + log2_pseudocount)` exactly once.
4. TCA exports one coordinate-preserving BED on the resulting log2 scale for
   each retained major cell type.
5. A scatter-input task validates the BED inventory against the WDL
   `Array[File]`. It emits aligned arrays of display names, safe slugs, and BED
   files.
6. The workflow scatters by aligned array index. Each call to
   `eQTLPrepareData` receives one TCA BED through `Log2CpmBed`, the shared
   sample list, and the required additional covariates.
7. The expression workflow rank-normalizes one branch and scales the other
   branch directly. It does not apply `log2(CPM + 1)` to the TCA output.
8. Each cell type and each transform performs its existing independent WGCNA
   connectivity filtering. Sample sets can therefore differ between cell
   types and between INT and scaled branches.
9. Each branch calculates selected and complete phenotype PCs and merges the
   selected PCs with the supplied additional covariates.
10. A manifest task validates the collected arrays and writes one row per
    retained cell type.

If any scatter call fails, the complete workflow fails. The workflow will not
publish a partial manifest.

## Scatter identity and filenames

The deconvolution inventory is authoritative for cell-type display names and
slugs. Before the scatter, validation will require:

- One inventory row per WDL BED file.
- Exact agreement between inventory basenames and localized BED basenames.
- Unique, nonempty cell-type names.
- Unique, nonempty safe slugs.
- The expected `log2_cpm` scale.
- Positive gene and sample counts.

Each scatter call will use this prefix:

```text
<OutputPrefix>.<cell_type_slug>
```

For example:

```text
cohort.cd4_t_cells.expression.INT.bed.gz
cohort.cd4_t_cells.expression.scaled.bed.gz
```

The manifest and all top-level file arrays will retain deconvolution inventory
order.

## Manifest contract

The end-to-end workflow will write a tab-separated manifest with these columns:

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
| `int_connectivity_outliers` | Connectivity outlier report for the INT branch. |
| `scaled_connectivity_outliers` | Connectivity outlier report for the scaled branch. |

The TSV will contain stable basenames, not task-local paths. WDL tasks cannot
discover the final Terra or Cromwell output URI that is assigned after output
delocalization. The workflow will expose one ordered `Array[File]` for every
file-valued manifest column. These arrays are the authoritative WDL file
outputs.

The workflow will not add raw or residualized BED files to this manifest. The
scattered subworkflows can create their normal raw outputs internally, but the
end-to-end public contract focuses on the requested INT and scaled QTL-ready
files.

## Upstream provenance outputs

The end-to-end workflow will also expose the important standalone
deconvolution outputs:

- Original and combined proportions.
- TCA weights and fitted model.
- Cell-group filter report.
- Original cell-type BED array and inventory.
- TCA reconstruction metrics and QC plots.
- Effective parameters, logs, and deconvolution output manifest.

This keeps the scientific model and retention decisions auditable from the
same workflow submission.

## Runtime and container design

The migrated deconvolution tasks will use a dedicated image built from:

```text
envs/CellTypeSpecificExpression/Dockerfile
```

The Dockerfile will use a pinned Micromamba base. The environment will pin R,
dtangle, TCA, tidyverse, Bioconductor, testing, and other direct package
dependencies. It will copy the migrated scripts to one stable path in the
image.

The scattered `eQTLPrepareData` calls will continue to use the standard
`prepare_QTL` image. The workflow will use distinct image inputs for
deconvolution and eQTL preparation so that the specialized TCA dependency set
does not enlarge or destabilize unrelated QTL workflows.

The end-to-end workflow will expose one global preemptible-attempt value and
one global retry value. New and updated task interfaces will pass these values
through the deconvolution and cell-type eQTL path. Each task will retain its
task-specific CPU, memory, and disk settings.

All new WDL command blocks will log stage start, important input dimensions or
settings, output paths, and completion. Any existing WDL task that must change
to support the new scatter interface will receive equivalent command logging.

## Validation and failure behavior

The workflow will stop with a direct error when:

- The deconvolution BED inventory and BED array do not match.
- An inventory cell type, slug, or path is empty or duplicated.
- The inventory scale is not `log2_cpm`.
- The pseudocount is negative or non-finite.
- A zero CPM value is present while the pseudocount is zero.
- A scatter output array has the wrong length or order.
- A required INT, scaled, PC, covariate, or outlier file is absent.
- The manifest would contain a task-local directory instead of a stable
  basename.
- The required additional covariates cannot be merged for a retained cell
  type.

Each scatter call owns its own output prefix. This prevents filename collisions
when calls run concurrently.

## Dockstore and documentation

`.dockstore.yml` will register both public entry points:

- Standalone cell-type deconvolution.
- End-to-end cell-type-specific eQTL preparation.

The workflow catalog and user documentation will describe:

- Linear-CPM input requirements.
- The one-time `log2(CPM + log2_pseudocount)` transformation in
  deconvolution.
- The zero-default pseudocount and its effect on zero-value validation.
- The lack of a second log transformation in QTL preparation.
- Independent connectivity filtering.
- Manifest fields and authoritative WDL arrays.
- Terra input examples for LM22 and precomputed-proportion modes.

After the integrated workflow passes GitHub Actions, the old
`CellTypeDeconvolution` repository will receive a final README deprecation
notice that directs users to the new canonical paths in `prepare_QTL`. The old
repository will not be deleted or kept synchronized.

## Testing strategy

Local tests will cover logic that does not require the specialized container:

- Inventory and BED-array alignment.
- Cell-type name and slug validation.
- Scatter order and output-prefix construction.
- Manifest schema, stable basenames, and array alignment.
- Missing, duplicate, malformed, or out-of-order records.
- WDL 1.0 parsing and import resolution.
- Dockstore descriptor registration.
- No second log transform for `Log2CpmBed`.
- Pseudocount validation, zero-value behavior, and identical dtangle/TCA
  pseudocount use.
- Required WDL logging fields.

GitHub Actions will build the dedicated Micromamba image and run small
end-to-end smoke workflows. CI will cover both proportion modes and verify:

- dtangle and TCA dependency versions.
- One QTL-preparation scatter call per retained cell type.
- Independent INT and scaled outputs.
- Phenotype-PC and merged-covariate files.
- Connectivity-outlier reports.
- Exact manifest rows, columns, order, and basenames.
- Important upstream deconvolution outputs.

No local Docker build is required. The current local workstation lacks
`biomaRt`, so the existing expression integration test cannot run locally at
the design baseline. The WDL tests, phenotype-PC tests, methylation integration
test, and Dockstore validation pass. GitHub Actions will provide the full
container-dependent verification.

## Branch and merge strategy

Development will occur on `codex/cell-type-specific-expression`, based on
`feat/log2-cpm-bed-input`. This creates a temporary stacked dependency on the
pre-normalized-expression work.

The intended final state is that both branches merge into `main`:

1. Merge the focused pre-normalized-expression branch.
2. Rebase or retarget the cell-type-specific branch onto `main` without
   duplicating the already merged commits.
3. Merge the cell-type-specific integration after its checks and review pass.

The worktree will remain available during pull-request review.
