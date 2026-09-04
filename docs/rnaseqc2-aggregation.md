# Batched RNA-SeQC aggregation

This workflow merges RNA-SeQC files in two stages. Each scatter shard merges at most 100 samples by default. The final task reads one row from each batch GCT and writes one cohort row. It does not load a complete GCT into memory.

## Files

- [`workflows/expression/rnaseqc2_aggregate_batched.wdl`](../workflows/expression/rnaseqc2_aggregate_batched.wdl): WDL 1.0 workflow for Terra.
- [`scripts/expression/merge_rnaseqc.py`](../scripts/expression/merge_rnaseqc.py): standard-library Python merger.
- [`envs/RNASeQCAggregation/`](../envs/RNASeQCAggregation): dedicated image and pinned package versions.
- [`sample_manifest.example.tsv`](../workflows/expression/examples/rnaseqc2_aggregation/sample_manifest.example.tsv): combined manifest example.
- [`inputs.example.json`](../workflows/expression/examples/rnaseqc2_aggregation/inputs.example.json): example Terra input JSON.
- [`tests/rnaseqc2_aggregation/`](../tests/rnaseqc2_aggregation): local unit and static WDL tests.

## Input manifests

The workflow uses one tab-separated sample manifest. Do not quote fields or add empty lines. Use one of these exact headers. In the examples below, `\t` means a tab character; the example TSV file contains actual tabs.

Without insert-size files:

```text
sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv
```

With insert-size files:

```text
sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\tinsert_size_hist
```

Each data row identifies one sample and all of its RNA-SeQC files:

```text
GTEX-001\tgs://bucket/sample-1/shared.gct.gz\tgs://bucket/sample-1/shared.gct.gz\tgs://bucket/sample-1/shared.gct.gz\tgs://bucket/sample-1/shared.tsv.gz
```

The path columns are:

- TPM GCT
- gene-read GCT
- exon-read GCT
- RNA-SeQC metrics TSV
- optional insert-size histogram

Each path must be a `gs://` URI. Repeated source basenames are allowed. The workflow stages every file with its manifest row number, sample ID, and file type. For example, it can stage two files named `shared.gct.gz` as `000001.GTEX-001.gene_tpm.gct.gz` and `000002.GTEX-002.gene_tpm.gct.gz`.

Each URI must include a bucket and object path. Keep the referenced objects immutable during a run. Terra fingerprints the manifest file for call caching, but it does not fingerprint objects that are listed as text inside the manifest.

Sample IDs must be unique. They can contain letters, numbers, periods, underscores, and hyphens. The workflow uses these IDs in all aggregate output tables, even when the individual RNA-SeQC files contain a different internal sample name.

The insert-size column is optional. If you include the column, every data row must contain an insert-size URI. If you omit the column, the workflow does not create an insert-size output.

## Terra use

Upload the combined sample manifest to a bucket that the Terra workspace can read.
The image includes `merge_rnaseqc.py`; do not upload a separate script.
When upgrading from the previous workflow, remove the `merge_script` input from
the Terra configuration or input JSON.

Import `rnaseqc2_aggregate_batched.wdl` into Terra. Set `batch_disk_space_gb` for the largest input batch, its temporary files, and its four or five batch outputs. Set `merge_disk_space_gb` from the total compressed size of all batch outputs and the expected final outputs. Include extra space for localization and temporary output files.

The output `prefix` must start with a letter or number. It can also contain periods, underscores, and hyphens.

The workflow uses `ghcr.io/aou-multiomics-analysis/prepare_qtl-rnaseqc2-aggregation`.
The default `docker_image` is pinned to the tested image's SHA-256 digest.
Updating the registry's `main` tag does not change this workflow version.
All three tasks use the same image and call the bundled script at
`/opt/prepare_qtl/scripts/expression/merge_rnaseqc.py`. If you override
`docker_image`, use a tested image digest that contains this script.

The image uses micromamba with pinned Python, `gsutil`, `crcmod`, Bash,
coreutils, findutils (`xargs`), and gawk packages from conda-forge. It does not
include RNA-SeQC, STAR, RSEM, R, or scientific Python libraries. It merges
existing RNA-SeQC outputs; it does not generate per-sample QC files.

The image contains no cloud credentials. `/etc/boto.cfg` tells standalone
`gsutil` to obtain credentials from the service account attached to the Terra
task VM. The account needs read access to all objects in the
manifest. This has not been tested with a live Terra submission.

The insert-size output is an `Array[File]` with zero or one file. WDL 1.0 uses this form so that the output can be absent without creating an invalid empty gzip file.

The metrics output is a cohort table. It has one row per sample and a first column named `sample_id`.

## Validation

The merger stops when it detects any of these conditions:

- duplicate output sample IDs;
- an incorrect manifest header or column count;
- an invalid sample ID or non-`gs://` path;
- different GCT dimensions;
- a different gene or exon at the same row;
- different sample order between GCT and metrics tables;
- different metric names or metric order;
- malformed insert-size bins or unsorted batch-level insert-size bins;
- an incomplete or extra GCT row.

The merger writes to a temporary file and publishes the result only after validation succeeds.

## Local checks

These checks do not build or run a Docker image:

```bash
python3 -m unittest discover -s tests/rnaseqc2_aggregation -v
miniwdl check workflows/expression/rnaseqc2_aggregate_batched.wdl
```

## Container checks and publishing

[GitHub Actions](../.github/workflows/rnaseqc2-aggregation.yml) builds the
Linux AMD64 image. No local Docker build is needed. The checks run the merger
tests against the script inside the image, with only the tests mounted from
the repository. They also check compiled checksum support and local file
copying with `gsutil`. A credential-selection check replaces only the metadata
server call and verifies that `gsutil` creates VM service-account credentials.
Real WDL validation and final-merge tasks run with
synthetic files, with and without insert-size tables. These checks do not
access private buckets or submit jobs to Terra.

Only a trusted branch push publishes the tested image to GHCR, with a
`sha-<full-commit>` tag. A push to `main` also updates the `main` tag.
Pull-request and manual runs build and test without publishing. Each run
stores the resolved package list as an artifact. The GHCR package must permit
anonymous pulls before Terra can use it without registry credentials.
