# CIBERSORTx fraction estimation on Terra

This WDL 1.0 workflow estimates relative cell fractions using an existing signature matrix. It supports optional S-mode correction and produces compatible input files for the separate HiRes workflow. It does not create a new signature matrix from single-cell data.

The complete Fractions workflow has **not been tested on Terra**. No cloud job was submitted. The user reported a successful Terra test of the separate HiRes workflow; that does not validate this workflow.

## Inputs

Use unadjusted linear expression and a compatible unadjusted signature matrix. For RNA-seq, leave `quantile_normalization=false`. Do not supply the previously adjusted signature with an unadjusted mixture by default.

| Input | Meaning | Default |
|---|---|---|
| `mixture` | Gene-by-sample expression matrix; first column contains gene symbols | Required File |
| `signature` | Gene-by-cell-type signature; first column contains gene symbols | Required File |
| `username`, `token` | CIBERSORTx account credentials | Required Strings |
| `smode` | Apply S-mode correction before fraction estimation | `false` |
| `refsample` | Reference profiles with replicates | Optional File; required here for S-mode |
| `source_geps` | Source expression profiles with the signature's cell labels and order | Optional File; required here for S-mode |
| `quantile_normalization` | Apply quantile normalization | `false` |
| `permutations` | Permutations for the native fit P-value | `100` |
| `threads` | Requested CPUs and thread limit for numerical libraries | `8` |
| `memory_gb`, `disk_gb` | Terra task memory and disk requests | `16`, `30` |

For a first operational test, set `permutations=0` to omit permutation testing. This does not establish meaningful fit P-values. The native Fractions program has no `--threads` flag. The workflow sets thread limits in the environment but cannot guarantee parallel fitting or a speed improvement from additional CPUs.

The native program can default `sourceGEPs` to the signature. This wrapper requires explicit `refsample` and `source_geps` for S-mode to use the project’s known reference triplet. The optional flags are passed only when the corresponding File input is present.

## Run in Terra

1. Add `workflows/cell_type_specific_expression/cibersortx_fractions.wdl` as a workflow method.
2. Upload the source files to a workspace-accessible bucket.
3. Replace the `gs://YOUR_BUCKET` and credential placeholders in the example JSON. These JSON files are workflow input examples, not script argument wrappers.
4. Import the inputs and submit the workflow in Terra on an x86 VM.

[`inputs.example.json`](../examples/cibersortx_fractions/inputs.example.json) uses the existing 100-sample unadjusted mixture and signature without correction. [`inputs.smode.example.json`](../examples/cibersortx_fractions/inputs.smode.example.json) uses the existing 50,000-cell Tabula reference triplet for correction:

- `mixture_100_samples_linear_CPM.txt`
- `signature_matrix_tabula_sapiens_50000.txt`
- `refsample_tabula_sapiens_50000.txt`
- `sourceGEP_tabula_sapiens_50000.txt`

These example filenames describe the project reference files; the matrices are user-supplied and are not included in this repository. For the 500-sample experiment, replace the mixture with `mixture_500_samples_linear_CPM.txt`. Fractions will be estimated for those 500 samples; the HSPE fractions are not an input to this workflow.

## Outputs for HiRes

Map these Fractions outputs to the existing HiRes WDL:

| Fractions output | HiRes input |
|---|---|
| `mixture_for_hires` | `mixture` |
| `signature_for_hires` | `signature` |
| `fractions_only` | `fractions` |

Supply the gene-subset File separately to HiRes. `signature_for_hires` retains only genes present in the corresponding mixture. This avoids the absent-signature-gene error from the earlier HiRes test. Fraction-only rows follow input sample order and columns follow signature cell order. Values are preserved and diagnostic columns are removed. Relative fractions must be finite, within [0,1], and sum to one within 0.0001.

`fractions` preserves the native table, including P-value, Correlation, and RMSE when present. `adjusted_mixture` and `adjusted_signature` are arrays with one File for S-mode and no Files otherwise. S-mode negative adjusted expression values are preserved and counted. `input_validation`, `output_validation`, `run_log`, `task_stdout`, and `task_stderr` provide checks and logs.

On failure, use the Terra task stdout and stderr links. Successful output fields can be absent when the task fails.

## Localization and logging

All input data files are typed as File or File? through workflow and task scope. The task verifies localized readability and rejects unresolved cloud URIs before native computation. It stages files in `/src/data`, sets a writable `/src/outdir`, and calls `/src/CIBERSORTxFractions` from `/src`.

Strings are passed as quoted named arguments, then as a native process argument list. No shell evaluation is used for the native command. The log reports start and completion times, matrix dimensions, validation, native activity, exit status, output checks, and output filenames. Credential option lines and literal credential text are omitted from the native log. Credentials still appear in workflow inputs and the rendered command; use the workspace's access controls for those artifacts.

The official amd64 image is pinned to:

```
cibersortx/fractions@sha256:9dc06b0a3f58d12a81cc962c9d2147b2b5edb6743f44dc2ac6d3f59fe7418edc
```

The digest, architecture, and entry point were verified from cached image metadata. No new Dockerfile or image build is needed.

## Checks

Local checks use a pinned WDL checker and a small native-process fixture. They exercise the rendered command, simulated cloud-to-local File handling, optional inputs, unusual path and credential characters, credential filtering, native failures, output labels, fractions, and shared signature genes. These checks do not run the authenticated CIBERSORTx algorithm.

```
python -m WDL check workflows/cell_type_specific_expression/cibersortx_fractions.wdl
python scripts/check_wdl_file_scope.py workflows/cell_type_specific_expression/cibersortx_fractions.wdl
python -m unittest discover -s tests/cibersortx_fractions -v
```

Use a Python environment with `miniwdl==1.14.2`. The static regression check rejects workflow-scope file writers in declarations, calls, outputs, scatters, and conditionals.

`.github/workflows/cibersortx-fractions.yml` runs these checks and a wrapper smoke test in the pinned image. The smoke test verifies the image's Python runtime, fixed paths, and output collection with a synthetic executable in place of the authenticated calculation. It requires no user data or token and does not build an image. The GitHub smoke test has passed for both plain mode and S-mode. It is not a complete Terra test or a test of scientific accuracy.
