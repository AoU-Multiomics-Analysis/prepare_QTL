# CIBERSORTx marker derivation on Terra

The WDL 1.0 workflow `CIBERSORTxMarkers` derives a marker signature from a labeled single-cell reference. It uses the marker derivation functions in the existing CIBERSORTx Fractions image. The default marker settings match the latest local Tabula Sapiens run with 50,000 cells: five replicates, sampling 1.0, fraction 0.0, 300 to 500 genes per cell type, and a q-value threshold of 0.01. Set `sampling` to `0.5` to use the earlier run's sampling value.

The complete workflow has **not run on Terra**. No cloud job or new authenticated marker derivation has been run for this addition. The GitHub Actions command test uses synthetic outputs and does not test the CIBERSORTx algorithm or biological results.

## Terra setup

1. Register `workflows/cell_type_specific_expression/cibersortx_markers.wdl` in a workflow method repository that Terra can use.
2. Upload the single-cell reference as a plain `.tsv` or gzip-compressed `.tsv.gz` file to a Google Cloud Storage bucket that the workspace can read.
3. Replace the bucket, username, and token placeholders in `examples/cibersortx_markers/terra.inputs.json`. Import this file into the Terra workflow configuration.
4. Review the resource settings and submit the workflow from Terra. Use an x86 VM for the amd64 image.

The example uses `single_cell_reference_all_cells_relabelled_no_globins_bulk_detected.tsv.gz`. This prepared input contains 22,278 genes and all 70,458 accepted cells from the reviewed Tabula Sapiens pool. It preserves the reviewed CD8/NK labels and the Job5 gene filters. It is stored in `analysis/cibersortx_relabelled_input_20260918/all_cells/` in the analysis workspace. These data files are supplied separately; they are not included in this repository. This input differs from the earlier 50,000-cell export, which used the original atlas labels. The workflow settings above come from that earlier command; the complete Job5 website settings were not recovered.

## Reference format

The `reference` input is a tab-delimited gene-by-cell matrix. The first column contains unique gene identifiers. Each remaining column contains one cell. Its header contains the cell-type label, so repeated cell-type column labels are expected. For example:

```text
GeneSymbol	B	B	B	CD4_T	CD4_T	CD4_T
GENE_A	12	14	11	1	0	1
GENE_B	0	1	1	8	9	10
```

Use finite, nonnegative expression values on the linear scale. Supply at least three cells for each cell type. Supply the full prepared single-cell reference, with consistent cell labels and gene identifiers. Do not use an existing signature matrix or a log-transformed matrix as the reference. The task checks the reference before it starts the native program.

The reference remains a WDL `File` until command rendering. Cromwell localizes it before the task checks and opens it. An unresolved cloud URI or unreadable local file causes a localization error. The example JSON supplies workflow inputs; it is not passed to the native program as an argument file.

The task detects gzip from the file contents. It decompresses gzip input on task disk before it checks the matrix and starts CIBERSORTx. The supplied file remains unchanged. Invalid or truncated gzip files cause an input error. Plain TSV input continues to work. No compression flag is needed. The input report records the supplied and uncompressed sizes and SHA-256 checksums. The disk request must allow space for both files and the native outputs; the compressed size alone is not the disk requirement.

## Settings

| Input | Default | Use |
|---|---|---|
| `reference` | Required `File` | Labeled single-cell matrix, plain TSV or gzip-compressed TSV. |
| `username`, `token` | Required `String` values | CIBERSORTx account credentials. |
| `replicates` | `5` | Number passed to `--replicates`. |
| `sampling` | `1.0` | Cell sampling value passed to `--sampling`. |
| `fraction` | `0.0` | Native expression-filter setting passed to `--fraction`. |
| `min_genes`, `max_genes` | `300`, `500` | Marker search limits passed to `--G.min` and `--G.max`. The number of unique output genes can differ from these values. |
| `q_value` | `0.01` | Threshold passed to `--q.value`. |
| `max_condition_number` | `999` | Native setting passed to `--k.max`. |
| `cpu` | `4` | Requested VM CPUs. This does not set a native thread flag. |
| `memory_gb` | `32` | Requested task memory in GiB. |
| `disk_gb` | `100` | Requested local disk capacity in GB. |

The task sets `--single_cell TRUE` and runs marker derivation without a bulk mixture. It does not estimate fractions or apply S-mode correction. The `fraction` input is a native filtering parameter. It is not an output cell fraction.

Sampling 1.0 uses all cells in each reference group. The generated replicate profiles do not represent independent biological donors. Review the resulting signature before using it for an analysis.

The task uses this image:

```text
cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc
```

It calls `/src/CIBERSORTxFractions` directly and uses the image's fixed input and output paths. The output directory points into the task working directory so Cromwell can collect the files. No new Dockerfile or image build is required.

## Outputs and logs

| Output | Contents |
|---|---|
| `signature_matrix` | Derived gene-by-cell-type marker matrix. |
| `source_geps` | Native source gene-expression profiles. |
| `reference_sample` | Native inferred reference matrix, with `number of cell types × replicates` expression columns. |
| `phenotype_classes` | Native inferred phenotype-class matrix. |
| `signature_heatmaps` | Native signature heatmap PDFs, when generated. |
| `input_validation` | Reference dimensions and input checks. |
| `output_validation` | Output dimensions and matrix checks. |
| `run_log` | Task messages and filtered native messages. |
| `task_stdout`, `task_stderr` | Task output streams. |

Keep the signature and reference outputs together. They describe the same derivation and can support a later CIBERSORTx fractions or S-mode workflow. This workflow returns native expression values; it does not assign new units to them.

The log records input checks, dimensions, start time, native exit status, output checks, and completion time. Credential option lines are removed from the native log. The username and token remain in the workflow inputs and rendered task command, which Terra can retain in execution metadata. Use the existing workspace access controls for these values.

For a failed task, use Terra's task log links to inspect stdout, stderr, and the execution directory. Successful output fields can be absent when a task fails.

## Validation

Run the syntax and local command checks in an environment with `miniwdl==1.15.0`:

```bash
miniwdl check workflows/cell_type_specific_expression/cibersortx_markers.wdl
python3 scripts/check_wdl_file_scope.py workflows/cell_type_specific_expression/cibersortx_markers.wdl
python3 scripts/check_wdl_logging.py workflows/cell_type_specific_expression/cibersortx_markers.wdl
python3 tests/cibersortx_markers/test_markers.py -v
```

The file-scope check rejects file-writing functions at workflow scope. The task tests simulate cloud File localization and execute the rendered command with a synthetic replacement for the native calculation. They test file handling, command arguments, failures, and output checks.

The GitHub Actions job repeats the command tests inside the pinned Fractions image on an x86 runner. It pulls the image and does not build an image. It needs no CIBERSORTx credentials. Its synthetic files must not be used as analysis results. Syntax checks, simulated localization, and container command checks do not establish success on Terra. A complete Terra run remains untested.
