# CIBERSORTx HiRes on Terra

This WDL 1.0 workflow wraps the existing CIBERSORTx HiRes image. It is prepared for the 100-sample ZNF804A test with nine cell types and 11 genes. The default CPU request and `--threads` flag are both 8. The sampling settings are both 1. These low sampling settings are for a first test. They do not establish stable expression estimates.

The complete workflow has **not run on Terra**. No cloud job has been submitted. The local HiRes runs were left unchanged.

## Run on Terra

1. Upload `workflows/cell_type_specific_expression/cibersortx_hires.wdl` to a workflow method repository that Terra can use.
2. Upload the four files below to a Google Cloud Storage bucket that the workspace can read. These data files are supplied separately; they are not included in this repository.
3. Replace the bucket and credential placeholders in `examples/cibersortx_hires/terra.inputs.json`. Import the values into the Terra workflow configuration.
4. Review the inputs and submit the test in Terra. Use an x86 VM for the amd64 image.

| WDL input | Source file |
|---|---|
| `mixture` | `mixture_100_samples_Smode_adjusted.txt` |
| `signature` | `signature_matrix_Smode_adjusted_shared_genes.txt` |
| `fractions` | `fractions_100_samples_Smode_fraction_only.txt` |
| `gene_subset` | `gene_subset_ZNF804A_controls.txt` |

Username and token are typed String inputs. The command passes them to `--username` and `--token`. The four data inputs remain typed File values until command rendering. Cromwell must localize them before the task starts. The example JSON is a workflow input example. It is not used as an argument wrapper for a script.

The signature file contains 3,023 shared genes. The 136 absent signature genes were removed for the local test. The fraction file contains only the nine cell fractions. Do not supply the original file with P-value, Correlation, and RMSE columns. Those columns caused the earlier matrix dimension error.

The fractions file can contain more samples than the mixture expression file. The task selects fraction rows by sample ID and puts them in the same order as the mixture sample columns. It preserves the fraction values and leaves the input file unchanged. Every mixture sample must occur once in the fractions file. Missing samples and duplicate labels cause the task to stop before HiRes starts. The selected rows must pass the fraction checks, including a positive total for each cell type.

The inputs already have S-mode adjustment. The wrapper uses `--QN FALSE` and does not repeat batch correction. It preserves negative values in the adjusted matrices and counts them in the input report. It does not add absent genes or change expression values.

## CPUs and image paths

Change `threads` to request more CPUs. For example, 16 requests 16 CPUs and passes `--threads 16`. More threads may reduce run time. Speed depends on the work that HiRes can run in parallel. Increasing threads beyond the work available for an 11-gene test may provide little benefit.

The wrapper uses the existing image pinned to this digest:

```text
cibersortx/hires@sha256:e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e
```

The image has an amd64 binary and an entry point of `./CIBERSORTxHiRes`. The task command calls that binary directly from `/src`. It stages the localized inputs under `/src/data`. It writes the selected fraction rows to `/src/outdir/fractions.txt`, because HiRes reads `--cibresults` there. Both directories point to the task's working directory. Cromwell can collect the outputs from that directory. No new image or Dockerfile is required. See [Cromwell container commands](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/) and [CPU and resource attributes](https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/).

## Logs and outputs

The command prints input checks, the run start, the exit status, and output checks. It prints a status line every minute while the native command is active. This line reports elapsed time. It does not prove that the fit is making progress or give a finish time. Native messages are included. Account and token option lines are omitted from the log.

Successful task outputs are:

- `expression_matrices`: one sample-level matrix per cell type.
- `input_validation`: sample IDs, cell types, gene counts, sampling settings, negative-value counts, and original input file checksums. The fields `fraction_input_sample_count`, `fraction_retained_sample_count`, and `fraction_dropped_sample_count` record the fraction row selection. These counts also appear in the run log.
- `output_validation`: sample and gene checks, missing-value counts, and counts of values equal to 1.
- `run_log`: `hires.log`, with task and native messages.
- `task_stdout` and `task_stderr`: the task streams.

For this input set, the output check requires nine matrices, the same 100 sample IDs in the same order, and all 11 subset genes. Missing estimates are allowed and counted. The supplied HiRes README documents 1 and NA as estimates with insufficient expression evidence or statistical power. The wrapper preserves them. Treat the native expression scale as unverified. It is not labeled CPM.

For a failed task, use Terra's task log links to read stdout, stderr, and the execution directory. Successful output fields may be absent when the command fails.

## Checks performed

`miniwdl check` passed. The tests execute the rendered task command. They simulate cloud File localization, including paths with spaces and apostrophes. They check literal credential flags, unresolved cloud paths, invalid fraction columns, output sample IDs, and native failure status. A static syntax-tree check rejects workflow-scope file-writing functions. These checks do not establish Terra integration.

The CIBERSORTx HiRes GitHub Actions check runs these tests locally and inside the pinned image on an x86 runner. It pulls the existing image and replaces the native calculation with a synthetic executable. It needs no account token and does not build an image. Changes to this WDL, its tests, its input example, or the file-scope check select this job. The complete native algorithm is not tested by this job.

The existing image also passed a task-only smoke test with the real four input files. A synthetic executable replaced the authenticated native calculation. The test checked 100 samples, nine cell types, 11 genes, fixed staging paths, CPU flags, logs, and output validation. Its matrices are synthetic and must not be used for analysis. This test did not run the CIBERSORTx algorithm.

To repeat the static checks, use a Python environment with `miniwdl==1.15.0`:

```bash
miniwdl check workflows/cell_type_specific_expression/cibersortx_hires.wdl
python3 scripts/check_wdl_file_scope.py workflows/cell_type_specific_expression/cibersortx_hires.wdl
python3 -m unittest discover -s tests/cibersortx_hires -v
```

After a real run completes, check gene expression variation and correlations across cell types. Exclude the documented insufficient-estimation markers from usable estimates. Do not interpret correlations between nearly constant estimates as evidence of biological sharing.
