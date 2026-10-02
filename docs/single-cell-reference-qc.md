# Tabula Sapiens reference filtering

This directory preserves the code used to remove predicted doublets, review
cell labels with SingleR, and select cells for the CIBERSORTx reference.
The scientific scripts are unchanged from the local analysis. These are
**cohort-specific source snapshots**, not a new general-purpose workflow.
Keep the input files outside the repository.

The code is in
[`workflows/cell_type_specific_expression/reference_qc`](../workflows/cell_type_specific_expression/reference_qc).
[`source_manifest.json`](../workflows/cell_type_specific_expression/reference_qc/source_manifest.json)
records each original source path and its SHA-256 checksum.

## Method and order

1. Start with the existing labeled Tabula Sapiens blood reference. Its gene
   set already has the study/Polaris bulk-detection filter, the 17 explicit
   gene exclusions including globins, and the X/Y gene exclusions. These
   scripts preserve those upstream rules; they do not reconstruct the
   upstream gene-selection or CD8/NK relabeling steps.
2. Export full raw UMI counts from the original H5AD for each eligible
   capture. Score all 71,101 cells in 22 captures from eight donors.
   The 643 context cells outside the 70,458-cell reference are scored but
   are not added to the reference. TSP10 remains excluded.
3. Run scDblFinder separately for each capture. Use `clusters = FALSE`,
   the automatic expected doublet rate, the package defaults, and the
   saved per-capture seeds. Remove reference cells called doublets.
4. Review the retained singlets with SingleR. Use Monaco main, Monaco
   fine, and HPCA main independently. Classify each capture, then prune
   scores once across all retained cells for each model. Sum duplicate
   gene symbols. Exclude X/Y-associated features for this review, but
   **keep globins in the raw counts used for label review**.
5. Filter the reference using the Monaco comparison table. Keep current
   CD8 cells only if their Monaco fine call supports CD8 and passed
   pruning. Remove plasma cells only if both Monaco models are uncertain.
   Preserve the other seven cell groups, including erythroid and platelets.

The matrix export preserves raw integer UMI counts, cell labels, and cell
order. It does not normalize counts or apply a log transform. It removes
gene rows that become zero across every retained cell and records them.
It writes gzip directly and verifies the decompressed output by checksum.

Monaco main and fine share one reference dataset. Their agreement is not
independent confirmation. A SingleR disagreement does not prove that a
cell is poor quality. The CD8 selection is restrictive and can remove
real CD8 states as well as incorrect labels.

## Source layout

| Directory | Purpose |
|---|---|
| `scdblfinder/` | Original local Snakefile, count export, R doublet detection, summaries, plots, and reference export |
| `singler/` | Original local Snakefile, R classification, label-coverage rules, and review summaries |
| `singler_filter/` | Final cell-selection and gzip matrix-export scripts and their tests |

The two Snakefiles are **local Snakemake workflows**. They are not WDL
workflows and have not run on Terra. No Dockerfile or image build is
included. CIBERSORTx marker derivation, S-mode correction, and fraction
estimation are separate steps and are not part of this code.

## Runtime

The completed analysis used R 4.3.1 and Bioconductor 3.17, with
scDblFinder 1.14.0, SingleR 2.2.0, celldex 1.10.1, and
SingleCellExperiment 1.22.0. The R scripts also need Matrix, readr,
dplyr, tibble, BiocParallel, and the packages loaded by those tools.
The doublet plot needs tidyverse and scales.

The scientific Python scripts used Python 3.14.6, numpy 2.5.3,
pandas 2.3.3, scipy 1.17.1, and h5py 3.16.0. The controller used
Snakemake 8.30.0 under Python 3.11.16. The executable and private R
library paths are supplied in the local configuration. Set `python_lib`
to an additional Python library directory if needed, or to an empty
string when the selected Python environment already has the packages.

[`runtime_versions.json`](../workflows/cell_type_specific_expression/reference_qc/runtime_versions.json)
records the observed direct versions. It is not a complete environment
lock. A newly installed environment or updated reference RDS file is
not guaranteed to reproduce all calls from the original analysis.

## Inputs

Replace every `/path/to/` value in the example configurations with a
readable absolute local path. The local JSON files configure Snakemake;
they are not WDL argument wrappers. Each script receives named CLI
arguments from its Snakefile.

| Input | Required content |
|---|---|
| Original H5AD | Tabula Sapiens Blood, with `obs`, `raw/var`, and raw CSR counts in `raw/X` |
| Selected-cell manifest | One row per reference column, with `cell_id`, `atlas_row`, `matrix_column`, `matrix_header`, `donor_id`, and `10X_run`; column numbers start at 2 |
| Capture plan | `capture`, `expected_cells`, and `seed`; all 22 captures must match the H5AD |
| Gene validation | Upstream sex-exclusion record, including source gzip checksum, gene counts, group counts, and prior gene exclusions |
| Reference matrix | Gzipped gene-by-cell TSV with raw nonnegative integer counts and repeated cell-type column labels |
| S100 review table | The original 226-cell review table; includes `cell_id` and `CD8_top10pct_S100` |
| Focused marker table | Earlier CD8 review with `cell_id`, `matrix_header`, `donor_id`, and the UMI/CPM marker fields used by `scdblfinder/scripts/summarize.py` |
| CLC review table | The six-cell marker-review table from the prior inspection; this is an external input, not an output of the copied doublet Snakefile |
| Capture exports | `counts.mtx.gz`, `features.tsv`, and `cells.tsv` in each capture directory from the doublet step |
| SingleR feature map | `feature_id`, `feature_symbol`, and Boolean `exclude_sex`; mapped to GENCODE v48, with any X/Y-associated locus excluded |
| SingleR references | Saved celldex Monaco and HPCA reference RDS files, using gene symbols and `logcounts` |

The 226-cell and six-cell review tables support the historical diagnostic
summaries. They are required by these snapshots even though they do not
define the general scDblFinder or SingleR algorithms. The source data,
cell-ID tables, reference RDS files, and matrices are not committed here.

### Cohort and checksum checks

These snapshots deliberately reject a different cohort. For example,
the count-export script requires 71,101 eligible raw cells and 70,458
reference cells; the SingleR script requires 65,191 retained singlets.
The final filter requires the verified original post-doublet gzip
checksum and expects 63,445 retained cells. These checks are part of
the code that produced the current reference.

**A fresh doublet rerun can produce a different gzip checksum even when
its decompressed counts match**, because the original exporter includes
gzip timestamps. The final filter therefore needs the original verified
post-doublet reference to run unchanged. To use a newly generated file,
first verify its counts and cell order, then update the explicit checksum
check in `singler_filter/scripts/prepare.py`. Do not bypass that check
without inspecting the new input. Adapting these scripts to another
cohort also requires reviewing the fixed-count assertions and summaries.

## Local commands

Run from the repository root. Create output directories before using
them as Snakemake working directories. Edit the example configurations
and save your bound copies outside the repository. Use the approved
execution process for the target where the analysis will run.

```bash
mkdir -p /path/to/doublets /path/to/label_review

snakemake \
  --snakefile workflows/cell_type_specific_expression/reference_qc/scdblfinder/Snakefile \
  --configfile /path/to/scdblfinder.config.json \
  --directory /path/to/doublets --cores 4 --printshellcmds

snakemake \
  --snakefile workflows/cell_type_specific_expression/reference_qc/singler/Snakefile \
  --configfile /path/to/singler.config.json \
  --directory /path/to/label_review --cores 4 --printshellcmds

python3 workflows/cell_type_specific_expression/reference_qc/singler_filter/scripts/prepare.py \
  --reference /path/to/verified_original_post_doublet_reference.tsv.gz \
  --selected /path/to/doublets/reference/selected_cells.tsv \
  --predictions /path/to/label_review/summary/all_cell_label_review.tsv.gz \
  --gene-map /path/to/gencode_v48_genes.tsv \
  --out-dir /path/to/final_reference
```

The gene map for the final filter uses `regulator_name_gencode` and
`regulator_chromosome`. It verifies that the already sex-excluded input
has no X/Y-associated symbols. The selected cells and predictions must
have identical cell IDs, labels, and row order.

## Outputs and observed run results

The doublet run writes capture scores and classes, method metadata,
session information, decision tables, and the post-doublet matrix and
cell manifests. In the original analysis it removed 5,267 reference
cells and retained 65,191 cells with 21,524 genes.

The SingleR run writes raw and pruned labels, score matrices, delta
scores, prediction RDS files, raw marker counts, comparison tables,
and session information. **Classification does not modify the matrix**;
the final filter is a separate script.

The final filter writes `filtered_reference.tsv.gz`, `selected_cells.tsv`,
`excluded_cells.tsv`, `counts_by_donor.tsv`, `input_gene_symbols.tsv`, and
`validation.json`. In the original analysis it retained 63,445 cells and
21,517 genes: 262 CD8 cells and 522 plasma cells. The other seven groups
were preserved. Independent full-matrix checks established that every
retained count matched the parent reference.

The delivered file was later renamed to
`tabula_sapiens_blood.polaris_filtered.doublets_removed.SingleR_filtered.autosomes.tsv.gz`.
This is a descriptive filename; the scripts still write the original
output name. Here `autosomes` refers to the X/Y exclusion used for that
file. It does not add a separate strict chromosome-1-to-22 filter.

## Verification

```bash
python3 tests/reference_qc/run_checks.py
Rscript -e 'paths <- list.files("workflows/cell_type_specific_expression/reference_qc", pattern = "[.]R$", recursive = TRUE, full.names = TRUE); invisible(lapply(paths, parse))'
```

The Python checks test reference-column alignment, preservation of raw
counts, reference coverage, CD8 selection, plasma uncertainty, and
preservation of the other groups. They also check source checksums and
example configuration consistency. GitHub Actions runs these checks
and R parsing without a local image build, external data, or credentials.
They do not rerun the scientific algorithms or establish Terra support.

The full scDblFinder and SingleR analyses completed locally before this
publication. The final cell-selection step also completed and was
verified. The subsequent CIBERSORTx marker calculation exhausted local
memory; it is not a failure of these filtering steps. A better CD8
signature or more accurate fraction estimates has not yet been shown.
