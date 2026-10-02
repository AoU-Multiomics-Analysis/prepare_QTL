# CIBERSORTx HiRes on Terra

This WDL 1.0 workflow splits a target-gene list, runs CIBERSORTx HiRes for each list, and merges one expression matrix per cell type. Every HiRes task receives the same full mixture, signature, fractions, and cohort. Only its `--subsetgenes` list changes.

The user reported a successful Terra test of the earlier workflow with one HiRes task. The chunked workflow has **not been tested on Terra**. Its split and merge checks do not establish that separate gene runs produce the same native estimates as one complete run. Sampling and gene-processing order can affect the estimates. Compare shared genes in a complete run and a chunked run before using chunked results for analysis.

## Run on Terra

1. Register `workflows/cell_type_specific_expression/cibersortx_hires.wdl` as a workflow method.
2. Upload the four input files to a bucket that the workspace can read. The data files are supplied separately; they are not included in this repository.
3. Replace the bucket and credential placeholders in `examples/cibersortx_hires/terra.inputs.json`.
4. Set `genes_per_chunk` and review the resources for each task. Submit the workflow on x86 VMs for the amd64 image.

The example uses the earlier 100-sample test with nine cell types and 11 genes:

| Input | Example file |
|---|---|
| `mixture` | `mixture_100_samples_Smode_adjusted.txt` |
| `signature` | `signature_matrix_Smode_adjusted_shared_genes.txt` |
| `fractions` | `fractions_100_samples_Smode_fraction_only.txt` |
| `gene_subset` | `gene_subset_ZNF804A_controls.txt` |

`gene_subset` is the full target list, with one unique gene per line and no header. For all-gene estimation, supply a list of all target genes in the mixture. Every target gene must occur in the mixture. Keep the full expression matrix as `mixture`; do not reduce it to one chunk's genes. Keep all cohort samples in each task. Splitting samples changes the sample windows and estimation context.

`genes_per_chunk` is the maximum number of genes in each subset list. Its default is `1000`. It must be a positive integer. The workflow creates `ceil(number of target genes / genes_per_chunk)` nonempty lists in the original order. The last list can be smaller. A value of `1` creates one list per gene. A value larger than the target count creates one list. For example, 2,501 target genes with `genes_per_chunk=1000` produce lists of 1,000, 1,000, and 501 genes. The 11-gene example produces one list with the default setting.

The four data inputs remain typed File values until command rendering. Cromwell must localize them before each task opens them. The example JSON supplies workflow inputs; it is not used as a script argument file. Username and token are String inputs passed to the native CLI.

The signature must contain only genes present in the mixture. The example signature has 3,023 shared genes; 136 absent genes were removed for the earlier test. The fractions file must contain only cell fractions, with cell columns in signature order. Remove P-value, Correlation, and RMSE columns before supplying it to HiRes.

The fractions file can contain more samples than the mixture. Each HiRes task selects rows by sample ID and puts them in mixture column order. It preserves the values and leaves the input file unchanged. Missing samples, duplicate labels, invalid fractions, or a cell type with no positive fractions stop the task before HiRes starts.

The example mixture and signature already have S-mode adjustment. The task uses `--QN FALSE` and does not repeat correction. It preserves and counts negative adjusted expression values. Do not pair these inputs with differently processed reference files without checking their compatibility.

## Resources and image

| Task | CPU | Memory | Disk |
|---|---:|---:|---:|
| `SplitGeneSubset` | 1 | 2 GiB | 1 GB |
| Each `RunHiRes` chunk | `threads`, default 8 | `memory_gb`, default 16 GiB | `disk_gb`, default 20 GB |
| `MergeHiRes` | 1 | `merge_memory_gb`, default 4 GiB | `merge_disk_gb`, default 20 GB |

The native settings and resource requests apply to **each chunk**, not to the complete scatter. Terra can run several chunks at once. Each task localizes the full input matrices, so more chunks can increase VM and storage costs. Set `memory_gb` for the full cohort and expression background, even when each target list is small. The split and merge tasks have separate resource requests.

The merger uses SQLite on task disk. Increase `merge_disk_gb` when the combined matrices are large. Allow space for the localized chunk matrices, the SQLite working copy, and the final matrices. Increase `merge_memory_gb` if needed for large headers or validation metadata. The merger streams expression rows instead of retaining the complete expression data in memory.

`threads` sets the CPU request and native `--threads` value for each HiRes task. The sampling settings default to `nsampling=1` and `nsampling2=1` for an operational test. These settings do not establish stable estimates. Use the same sampling settings for all chunks and for a complete-run comparison.

All tasks use the existing image pinned to:

```text
cibersortx/hires@sha256:e8da6850311d163e33a343d29a0d2ffc8b18c1ec4604995b8285c7bc2017c83e
```

Each native task calls `./CIBERSORTxHiRes` from `/src`. It stages localized inputs under `/src/data` and writes selected fractions to `/src/outdir/fractions.txt`, where HiRes reads `--cibresults`. Output directories point into the task working directory so Cromwell can collect the files. No new Dockerfile or image build is required.

## Merge checks and outputs

The merger checks exact target-gene coverage and rejects missing, repeated, unexpected, or overlapping genes. It also rejects different sample order, cell labels, input checksums, sampling settings, or window sizes across chunks. Matrices with identical names in different chunk directories remain separate File inputs. The merger restores the original target-list order in every final cell matrix.

The merger preserves expression text, including the native `1`, `NA`, `NaN`, and blank missing-value markers. HiRes documents `1` as insufficient evidence of expression or insufficient power for estimation. Do not treat these markers as measured expression. Native expression units remain unverified.

The existing workflow output names and types are retained:

| Output | Contents |
|---|---|
| `expression_matrices` | One merged sample-level expression matrix per cell type. |
| `input_validation` | Combined input checks and full target-gene context. |
| `output_validation` | Merged sample, gene, cell-type, and value checks. |
| `run_log` | Merge checks, native chunk logs, and completion messages. |
| `task_stdout`, `task_stderr` | Merge task streams. |

Additional outputs expose each stage:

| Output | Contents |
|---|---|
| `chunk_plan` | Target count, maximum list size, and chunk plan. |
| `chunk_gene_lists` | Generated subset lists, in chunk order. |
| `chunk_input_validations`, `chunk_output_validations` | Each native task's validation reports. |
| `chunk_logs` | Each native task's filtered log. |
| `chunk_stdout`, `chunk_stderr` | Each native task's streams. |
| `split_log` | Split task checks and completion messages. |

The logs report input checks, dimensions, start time, native exit status, output checks, and completion. Native tasks print an activity message every minute. This message is not a progress percentage. Credential option lines are omitted from native logs. Credentials remain in workflow inputs and rendered commands; use workspace access controls for those artifacts.

If a task fails, inspect its Terra stdout, stderr, and execution directory. Successful output fields can be absent after failure. Use the chunk reports and logs to identify a failed native call.

## Checks

Local task tests use synthetic native outputs. They check command arguments, full input forwarding, gene-list size, exact merge coverage, ordering, missing-value preservation, native failure status, and cloud-to-local File handling. The localization checks include generated chunk Files and task-local File lists. These tests do not execute the authenticated CIBERSORTx calculation.

Use an environment with `miniwdl==1.15.0`:

```bash
miniwdl check workflows/cell_type_specific_expression/cibersortx_hires.wdl
python scripts/check_wdl_file_scope.py workflows/cell_type_specific_expression/cibersortx_hires.wdl
python scripts/check_wdl_logging.py workflows/cell_type_specific_expression/cibersortx_hires.wdl
python -m unittest discover -s tests/cibersortx_hires -v
```

The static file check rejects file-writing functions at workflow scope. The logging check covers the split, native, and merge tasks. The Cromwell Womtool test includes this descriptor and checks imports and types; it does not run a Terra workflow.

`.github/workflows/cibersortx-hires.yml` runs all HiRes tests on the host and inside the pinned image on an x86 runner. It pulls the image and replaces the authenticated native calculation with a synthetic executable. It requires no user data or account token and does not build an image. Container checks do not establish native chunking equivalence or complete Terra compatibility.
