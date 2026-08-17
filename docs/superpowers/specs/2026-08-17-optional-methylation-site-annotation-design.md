# Optional methylation site annotation design

## Goal

Allow methylation cohort workflows to skip site annotation when it is not needed, while preserving current behavior by default and increasing annotation resources for enabled runs.

## Design

- Add `Boolean AnnotateSites = true` to the public methylation cohort workflow inputs.
- Conditionally call `Annotation.AnnotateMethylationCohortSites` only when `AnnotateSites` is true.
- Expose `PassingSiteAnnotations` as `File?` because no annotation file exists when annotation is disabled.
- Set the public annotation defaults to `AnnotationMemoryGB = 256` and `AnnotationDiskGB = 200`.
- Preserve all aggregation, connectivity, phenotype-PC, QTL-covariate, and other outputs regardless of annotation selection.

## Scope

Update the array-based cohort workflow, the manifest-based cohort wrapper, and the merge workflow so the option and resource defaults propagate consistently. Add WDL validation and regression checks for the default-enabled and disabled call structure.
