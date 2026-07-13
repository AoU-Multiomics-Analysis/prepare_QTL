# Workflow catalog

Workflows are grouped by data modality. Shared molecular-QTL helpers live in
`common/`; each other directory owns its public entry points and private
implementation modules.

## Directory layout

| Directory | Contents |
| --- | --- |
| `common/` | Phenotype-PC calculation, covariate merging, and phenotype residualization shared across modalities. |
| `expression/` | Expression-QTL preparation. |
| `proteomics/` | Olink normalization and protein-QTL preparation. |
| `splicing/` | Splicing-QTL preparation. |
| `methylation/` | Per-sample processing, cohort aggregation, connectivity analysis, annotation, and methylation-QTL outputs. |
| `genotype/` | VCF preparation, PLINK conversion, genotype PCs, allele frequencies, and dosage extraction. |

## Public entry points

| Area | WDL | Use |
| --- | --- | --- |
| Expression | `expression/prepare_eQTL.wdl` | Prepare expression QTL phenotypes. |
| Proteomics | `proteomics/normalize_pQTL.wdl`, `proteomics/prepare_pQTL.wdl` | Normalize Olink data and prepare protein QTL phenotypes. |
| Splicing | `splicing/prepare_sQTL.wdl` | Prepare splice QTL phenotypes. |
| Methylation | `methylation/ProcessMethylationSample.wdl` | Per-sample Terra-table processing. |
| Methylation | `methylation/AggregateMethylationCohort.wdl` | Cohort aggregation from one compact `CohortManifest`. |
| Methylation | `methylation/merge_methylation.wdl` | Legacy source-BED manifest/shard entry point. |
| Genotypes | `genotype/prepare_VCF.wdl`, `genotype/convertVCF2Plink.wdl`, `genotype/calculateGenotypePCs.wdl`, `genotype/calculateAF.wdl`, `genotype/calculateGenotypeDosage.wdl` | Genotype preparation and derived outputs. |

Every public entry point and the independently callable common helpers are
registered in [`.dockstore.yml`](../.dockstore.yml) and checked in CI.

## Internal methylation stages

- `methylation/AggregateMethylationCohortArrays.wdl`: array-based cohort implementation used only inside `methylation/merge_methylation.wdl`.
- `methylation/cohort_aggregation.wdl`: cohort-manifest expansion, sample validation, per-chromosome merging, cohort-wide aggregation, and annotation orchestration.
- `methylation/connectivity.wdl`: preliminary phenotype PCs, covariate-adjusted CpG correlation, and sample-connectivity filtering.
- `methylation/annotation.wdl`: annotation of passing methylation sites, called by the aggregation workflow.
- `methylation/qtl_covariates.wdl`: final phenotype PCs and optional TensorQTL covariate assembly.

`methylation/AggregateMethylationCohort.wdl` is intentionally a thin public
orchestrator over the methylation stage workflows. Keep its workflow-level
inputs and outputs stable for existing configurations; place new implementation
tasks in the stage module that owns their data.

## Maintenance rule

When adding, renaming, or removing a public workflow, update `.dockstore.yml`,
this catalog, and the corresponding user-facing documentation in the same
change. The workflow-validation action checks descriptor existence, uniqueness,
and WDL parsing on every pull request and push to `main`.
