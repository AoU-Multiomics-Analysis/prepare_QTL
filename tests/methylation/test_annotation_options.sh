#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "${repo_root}"

workflow_files=(
    workflows/methylation/cohort_aggregation.wdl
    workflows/methylation/AggregateMethylationCohortArrays.wdl
    workflows/methylation/AggregateMethylationCohort.wdl
    workflows/methylation/merge_methylation.wdl
)

for workflow in "${workflow_files[@]}"; do
    miniwdl check "${workflow}" >/dev/null
done

for workflow in \
    workflows/methylation/AggregateMethylationCohortArrays.wdl \
    workflows/methylation/AggregateMethylationCohort.wdl \
    workflows/methylation/merge_methylation.wdl; do
    grep -Fq 'Boolean AnnotateSites = true' "${workflow}"
    grep -Fq 'Int AnnotationMemoryGB = 256' "${workflow}"
    grep -Fq 'Int AnnotationDiskGB = 200' "${workflow}"
    grep -Fq 'File? PassingSiteAnnotations' "${workflow}"
done

for workflow in \
    workflows/methylation/cohort_aggregation.wdl \
    workflows/methylation/AggregateMethylationCohortArrays.wdl; do
    grep -Fq 'if (AnnotateSites)' "${workflow}"
    grep -Fq 'PassingSiteAnnotations' "${workflow}"
done

echo "Annotation option checks passed"
