# REGENIE files from QTL preparation

The preparation workflows now produce REGENIE files beside their BED and
TensorQTL files. Use the same genotype cohort and participant IDs for both.

| Entry point | REGENIE phenotype branches |
| --- | --- |
| `prepare_eQTL`, `prepare_pQTL`, `prepare_sQTL`, `prepare_mQTL` | INT, scaled, raw |
| `AggregateMethylationCohort`, `merge_methylation` | INT, raw |
| `prepare_cell_type_eQTL` | INT, scaled, raw for each cell type |
| `PhenotypePCQTL` | The selected phenotype PCs, with those PCs removed from the covariates |

Set `RegeniePhenotypesPerFile` to change the batch size. The default is **100
traits per file**. This is an export setting, not a measured performance limit.
Each batch contains all retained samples and at most that number of traits.
The converter holds one batch in memory. Large exports can still need substantial
disk space. The export task uses the preparation workflow's existing disk setting
(`disk_space`, `PcDiskGB`, `AggregateDiskGB`, or `disk_gb`). Increase that setting
to fit the uncompressed phenotype TSVs as well as the localized BED and inputs.

For each branch, the workflow adds these outputs. Replace `Int` with `Scaled`
or `Raw` as needed. The phenotype-PC workflow uses `Pc`.

| Output | Contents |
| --- | --- |
| `IntRegeniePhenotypes` | `Array[File]` of uncompressed phenotype TSV batches |
| `IntRegenieCovariates` | Optional covariate TSV |
| `IntRegenieTraitMap` | Safe REGENIE name, original trait ID, chromosome, start, end, batch filename |
| `IntRegenieBatches` | Batch filename, comma-separated phenotype names, trait count, sample count |
| `IntRegenieExcludedSamples` | BED sample IDs absent from the covariate table |
| `IntRegenieLog` | Export status and counts |

Cell type outputs have one outer array per cell type, in the order of the
existing `cell_types` and `cell_type_slugs` outputs. Thus, `IntRegeniePhenotypes`
is `Array[Array[File]]`. Each other export has one value per cell type. The
optional covariate output has type `Array[File?]`.

The phenotype and covariate TSVs have one participant per row:

```text
FID IID pheno_04b... pheno_709...
0 participant_1 1.23000 NA
0 participant_2 -2e-08 0.7
```

The actual files use tabs. `FID` is `0`; `IID` is the exact BED sample ID.
For a PSAM that has only `#IID`, the REGENIE pipeline adds `FID=0` in its task.
Numeric text is retained; outer padding in covariate cells is removed. Missing phenotypes (`NA`, `NaN`, `.`, or an empty cell)
become `NA`. The converter does not transform, impute, or residualize values.
Each original trait ID gets a stable name of the form `pheno_<24 hex digits>`.
Use the trait map to join association results to genes, proteins, splice events,
or CpG sites. Coordinates are copied from the input BED without conversion.

INT and scaled exports use the corresponding merged TensorQTL covariates when
`AdditionalCovariates` is supplied. Otherwise, they use the selected phenotype
PCs. Raw exports use only `AdditionalCovariates`; they do not use phenotype PCs
from a different transformation. If there are no covariate columns, the covariate
output is absent. Covariate names are kept and must be valid REGENIE names.

Samples follow BED order. With covariates, only samples present in both files
are exported. The exclusion report lists removed BED samples. Missing covariate
cells cause an error. Categorical covariate values are kept as text; select their
column names in the REGENIE `categorical_covariates` input.

Use a batch with the [REGENIE pipeline](https://github.com/AoU-Multiomics-Analysis/regenie):

1. Set `phenotypes` to one File from `IntRegeniePhenotypes`.
2. Find that filename in `IntRegenieBatches`. Copy its `phenotype_columns` value
   to the REGENIE input of the same name.
3. Set `covariates` to `IntRegenieCovariates` when that output exists.
4. Set `trait_type` to `quantitative`. Leave `inverse_normal_transform=false`
   when using an INT export, because those values are already transformed.
5. Run each phenotype batch separately with the required PLINK 2 inputs and
   common-variant panel for step 1. The preparation workflow does not run GWAS.

The batch table contains filename metadata. It is not a file-localization list.
For WDL dependencies, pass the typed phenotype File outputs directly. Optional
residualized BEDs remain available in their existing format. To export one, run
the standalone `workflows/common/ExportRegenie.wdl` helper with that BED. Choose
covariates with care to avoid a second adjustment for the same predictors.

The helper accepts an optional covariate File and `CovariateFormat="sample"`
(`ID`, `IID`, or `sample_id` followed by predictor columns) or `"tensor"`
(covariates as rows, samples in the header). It uses an existing pinned REGENIE
image. The Python source under `tools/regenie_export/` is also embedded in the
WDL so that no remote script download or image rebuild is needed. Tests check
that both copies match.

GitHub Actions checks formats, WDL syntax, command localization, and synthetic
two-step REGENIE scans for all four molecular QTL types. A complete workflow has
**not been tested on Terra**. No Terra job is submitted by these checks.
