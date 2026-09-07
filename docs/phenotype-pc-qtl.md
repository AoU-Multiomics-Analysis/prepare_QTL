# Genetic association scan of phenotype PCs

A two-step WDL 1.0 workflow for Terra/Cromwell:

1. Select existing phenotype PC rows from a tensorQTL covariate matrix. Write them as BED phenotypes and retain every other row as an adjustment covariate.
2. Run tensorQTL trans against PLINK 2 genotypes using the remaining covariates.

This does not calculate new PCs, select PCs for removal automatically, or modify the original covariate file. Inspect the resulting PC–variant associations and their biology before deciding which PCs to omit from a later trans-QTL analysis. The saved P-value cutoff is not an automatic exclusion rule.

## Inputs

Use `workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl` with its `workflows/phenotype_pc_qtl/tasks/tensorqtl_trans_compat.wdl` import. Keep their relative paths when uploading or registering the workflow. The example Terra inputs JSON is a workflow submission example, not a script argument wrapper. Replace its example bucket URIs and PC IDs.

The covariate TSV has one covariate per row, its ID in column 1, and sample IDs across the header:

```text
ID       sample1  sample2  sample3 ...
genPC1   0.1      0.2      -0.3    ...
exprPC1  -0.2     0.5      0.1     ...
exprPC2  0.6      -0.1     0.3     ...
```

Supply `phenotype_pc_ids` as the exact row names of **all phenotype PCs to exclude from adjustment and test**. Example: `["exprPC1", "exprPC2"]`. Names are not inferred. Genetic PCs, sex, cohort indicators, and any other unselected rows stay in the covariate file. If a phenotype PC is not selected, it remains an adjustment covariate. Sample order and numeric text are preserved; values are not normalized again.

Inputs must be numeric and finite, with unique sample and covariate IDs and no constant rows. Encode categorical covariates before submission. The workflow rejects constant covariates (including a user-supplied intercept), missing PC IDs, and insufficient residual degrees of freedom by count. It does not test full matrix rank.

All genotype files are explicit WDL `File` inputs. PVAR must be uncompressed text. PSAM must have an IID column; every phenotype sample must be present. Extra genotype samples are allowed because tensorQTL selects phenotype samples. The task stages the three localized files under a common prefix, so their original names and directories need not match.

## Placeholder coordinates

Every PC has default BED coordinates `chrY 0 1`. Coordinates are placeholders, not genomic positions of PCs. TensorQTL's sparse trans scan filters same-chromosome pairs within 5 Mb of the phenotype position. Thus autosomal and chrX tests are retained, but nearby chrY variants can be excluded. The task logs a warning when the PVAR contains the placeholder chromosome.

To avoid coordinate-based exclusion entirely, use a chromosome label absent from your genotypes, such as `PC_SCAN`. This uses tensorQTL's BED reader; a separate external BED validator that requires canonical chromosome names may reject this label. Chromosome labels should be consistent with the PVAR when using chrY. No same-chromosome-wide exclusion is added by this wrapper.

## Outputs and resource settings

- `phenotype_pc_bed`: selected PC rows with placeholder coordinates.
- `remaining_covariates`: unselected covariate rows.
- `pc_selection`: one row per input covariate, with its assigned role.
- `trans_qtl_pairs`: sparse Parquet table of saved PC–variant associations.
- `preparation_log` and `trans_log`: task messages and errors.

Defaults: MAF >= 0.05, saved nominal P < 1e-5, CPU execution with 4 threads, 64 GB RAM, and 100 GB working disk for the trans task. Increase memory and disk for large genotype inputs; these are not genome-wide resource guarantees. The wrapper uses the upstream command's genotype-loading behavior and does not scatter by chromosome. It does not compute PC-level genome-wide adjusted P-values. Sparse output means variants that fail the saved P cutoff are absent, not untested.

## Upstream compatibility

Source requested by the user: https://github.com/evin-padhi/tensorQTL_trans (redirects to AoU-Multiomics-Analysis/tensorQTL_trans).

Source task: https://github.com/AoU-Multiomics-Analysis/tensorQTL_trans/blob/ab5aa7befaae3c6d0b5f20178cb587f10b4da160/tensorQTL_trans.wdl

The original is a legacy draft WDL with implicit workflow inputs and no explicit workflow outputs. This package imports a local WDL 1.0 compatibility copy of its task rather than the unchanged legacy workflow. Changes are explicit input/output blocks, named and safely quoted arguments, localized PLINK staging, input checks, logs, a saved P-value cutoff, CPU runtime settings, and a fixed image digest. Sparse mode is used because the original declared output is a sparse pairs file. Its unused trans FDR option and hard-coded GPU zone were omitted.

The existing upstream `latest` image resolved on 2026-09-07 to:

```text
gcr.io/broad-cga-francois-gtex/tensorqtl@sha256:f6efb9e592eb32c46cb75070be2769b34381d60cbb2709d2885771324abfe32a
```

No new Dockerfile or image build is needed. `tensorqtl_docker` can override the image if required. The Python helpers use only the standard library and are embedded in the WDL task commands, so Terra does not need a separate script input. After editing a helper, run `python tools/phenotype_pc_qtl/sync_wdl.py`.

## Validation

From the repository root, with miniwdl 1.15.0 installed:

```bash
miniwdl check workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl
python -m unittest discover -s tests/phenotype_pc_qtl -v
java -jar /path/to/womtool-87.jar validate workflows/phenotype_pc_qtl/phenotype_pc_qtl.wdl
```

Local checks validate the complete WDL import graph and run rendered task commands. Localization tests remap synthetic `gs://` File values into different local directories, include shell metacharacters in filenames, and check that missed localization fails before computation. TensorQTL itself is stubbed in that boundary test. A static AST check rejects workflow-scope file-writing functions. The only `write_lines` call creates a task-local list of PC names; no File inputs are serialized before localization.

`.github/workflows/phenotype-pc-qtl.yml` checks this workflow. It runs validation and a separate synthetic numerical scan inside the pinned upstream image. The numerical test checks recovery of a strong simulated association and chrY versus absent-chromosome filtering. It pulls the existing image; it does not build an image.

**The complete workflow has not been run on Terra. The GitHub Actions smoke test runs on the pull request; inspect its result before use.** No cloud jobs have been submitted. WDL validation and local command tests do not establish Terra execution or image-runtime validation.
