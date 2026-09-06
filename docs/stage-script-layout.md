# Stage script layout

Scripts use the stage directories registered in `ci/image-stages.yml`.
The move does not change calculations, arguments, input scales, or output formats.

| Source directory | Release stage |
| --- | --- |
| `scripts/cell_type_specific_expression/estimation/` | `cell_estimation` |
| `scripts/cell_type_specific_expression/fit/` | `cell_fit` |
| `scripts/cell_type_specific_expression/export/` | `cell_export` |
| `scripts/cell_type_specific_expression/downstream/` | `cell_downstream` |
| `scripts/expression/prepare/` | `expression` |
| `scripts/expression/rnaseqc/` | `rnaseqc` |
| `scripts/methylation/` | `methylation` |
| `scripts/splicing/` | `splicing` |
| `scripts/proteomics/` | `proteomics` |
| `rust/methylation_filter/`, `rust/methylation_merge/` | `methylation_rust` |

`scripts/common/` contains shared tools, not one exclusive stage. The current
four files retain the common, expression, proteomics, splicing, and methylation
consumers. A new common file needs an explicit consumer review. Cell-type modules
stay in `scripts/cell_type_specific_expression/R/`; their shared dependencies
remain explicit. `bootstrap.R` stays at the cell-type root. Stage entrypoints
load it from their parent directory. Developer tools stay outside runtime stages.

New scripts inside a stage directory inherit its stage mapping. A new WDL task
must still receive the matching image input. CI checks its script paths and
image forwarding. This move does not create a separate Dockerfile per stage.

## Safe rollout

The release controller does not accept test or policy edits in an image-release
PR. Use two PRs; do not switch main's WDL commands to paths absent from its images.

1. **Preparation:** add the stage-directory copies and updated tests, CI filters,
   and maps. Retain the original entrypoints, Dockerfiles, WDL commands, and
   digest defaults. These original entrypoints are temporary compatibility
   copies, not a second maintained implementation. Do not use `release-ready`.
2. **Migration:** remove the original entrypoints and change the Docker copy
   paths and WDL commands. This PR must not change tests or release policy.
   Build and test the candidate images through the release workflow. Merge
   only after the new digest defaults are committed and checks pass.

Keep legacy source mappings during this release: Git reports a move as a
deleted path plus an added path. Both must have a stage. The `/tmp/PrepareExpression.R`
image entrypoint remains available; its canonical source is now in `prepare/`.
Other existing `/tmp/` entrypoints are unchanged.

No Terra run or local Docker build is part of this reorganization. Source tests
and descriptor checks do not establish that the full workflow has run on Terra.
