# Repository image routing

WDL tasks must not use a mutable tag for an image that this repository builds.
Each task runtime uses a declared `String` input. A workflow passes its image
input to every task or imported workflow that it calls.

## Image inputs

| Release family | Workflow input | Repository |
| --- | --- | --- |
| Expression | `DockerImage` | `ghcr.io/aou-multiomics-analysis/prepare_qtl` |
| Common | `DockerImage` | `ghcr.io/aou-multiomics-analysis/prepare_qtl` |
| Proteomics | `proteomics_docker_image` | `ghcr.io/aou-multiomics-analysis/prepare_qtl` |
| Splicing | `splicing_docker_image` | `ghcr.io/aou-multiomics-analysis/prepare_qtl` |
| Methylation | `methylation_docker_image` | `ghcr.io/aou-multiomics-analysis/prepare_qtl` |
| Methylation Rust | `methylation_rust_docker_image` | `ghcr.io/aou-multiomics-analysis/prepare_qtl-methylation-rust` |
| RNA-SeQC aggregation | `docker_image` | `ghcr.io/aou-multiomics-analysis/prepare_qtl-rnaseqc2-aggregation` |

The family-specific names keep release movement separate. For example, a
proteomics update does not change a methylation workflow default. When a family
workflow calls a common workflow or task, it passes its family input to the
common `DockerImage` input.

The methylation cohort workflow uses two image inputs. It passes the standard
image to R tasks and the Rust image to `MergeMethylationChromosome`. This rule
also applies through imported methylation workflows.

## Initial immutable references

The initial standard reference is the previously verified QTL reference:

`ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845`

The initial RNA-SeQC reference is:

`ghcr.io/aou-multiomics-analysis/prepare_qtl-rnaseqc2-aggregation@sha256:c2dc991dc99d8323fe6cc22375cd6560c131afc02122f8a6b2050eb7adba7652`

The initial methylation Rust reference is:

`ghcr.io/aou-multiomics-analysis/prepare_qtl-methylation-rust@sha256:16f631c34e0ce265d686335b91c18948607127178d95e7829070c97cd207d6ad`

On 2026-09-05, a read-only GHCR manifest request returned HTTP 200 and the
same `Docker-Content-Digest` for each exact digest reference. This check proves
that each manifest was published at that time. It does not test a task or a
complete workflow.

## Release-pin targets

`ci/release-pins.yml` maps each release family to the defaults that it can
change. A target without `scope` selects an input on the document workflow.
A task default uses this form:

```yaml
- path: workflows/common/ResidualizePhenotypes.wdl
  scope: task
  task: ResidualizePhenotypes
  input: DockerImage
```

The proposal parser must accept `scope: task` and require one exact `task`
name. It must select the named task input and reject an absent task, a duplicate
task name, a non-String input, or a non-literal default. If `scope` is absent,
the parser must keep the current workflow-input behavior. A target identity is
the tuple `(path, scope, task, input)`. The parser must reject duplicate target
identities.

The `common` release family owns shared task and workflow defaults. Family
entry points still pass their own image inputs into common calls. This design
prevents one shared default from replacing the independent expression,
proteomics, splicing, and methylation defaults.

## Validation boundary

`tests/test_repo_image_routing.py` loads the WDL documents with MiniWDL. It
checks task runtime inputs, immutable registered defaults, and evaluated call
input forwarding. Genotype WDL files are not in this policy because their
images are external and remain manually maintained.

These checks are local syntax and routing checks. The complete workflows have
not been tested on Terra. No cloud job was submitted for this change.
