# Workflow catalog

Public workflow descriptor paths remain flat in this directory so existing Terra method configurations, Dockstore registrations, and downstream imports remain valid.

## Public entry points

| Area | WDL | Use |
| --- | --- | --- |
| Expression | `prepare_eQTL.wdl` | Prepare expression QTL phenotypes. |
| Proteomics | `normalize_pQTL.wdl`, `prepare_pQTL.wdl` | Normalize Olink data and prepare protein QTL phenotypes. |
| Splicing | `prepare_sQTL.wdl` | Prepare splice QTL phenotypes. |
| Methylation | `ProcessMethylationSample.wdl` | Per-sample Terra-table processing. |
| Methylation | `AggregateMethylationCohort.wdl` | Cohort aggregation from one compact `CohortManifest`. |
| Methylation | `merge_methylation.wdl` | Legacy source-BED manifest/shard entry point. |
| Genotypes | `prepare_VCF.wdl`, `convertVCF2Plink.wdl`, `calculateGenotypePCs.wdl`, `calculateAF.wdl`, `calculateGenotypeDosage.wdl` | Genotype preparation and derived outputs. |

Every public entry point is registered in [`.dockstore.yml`](../.dockstore.yml) and checked in CI.

## Internal building blocks

These WDLs are imported by public workflows or support their implementation; they are intentionally not separate Dockstore entry points.

- `AggregateMethylationCohortArrays.wdl`: array-based cohort implementation used only inside `merge_methylation.wdl`; it avoids submitting large arrays through Terra's API.
- `calculate_phenotypePCs.wdl`, `MergeCovariates.wdl`, `ResidualizePhenotypes.wdl`: shared molecular-QTL helpers.

## Maintenance rule

When adding, renaming, or removing a public workflow, update `.dockstore.yml`, this catalog, and the corresponding user-facing documentation in the same change. The workflow-validation action checks descriptor existence, uniqueness, and WDL parsing on every pull request and push to `main`.
