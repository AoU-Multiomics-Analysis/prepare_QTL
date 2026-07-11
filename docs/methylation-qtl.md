# PacBio 5mC QTL workflow

[Back to main README](../README.md)

This guide describes how to prepare pb-CpG-tools site-level 5mC calls for molecular QTL mapping. The workflow applies coverage QC per sample, filters sites against the whole cohort, applies a cohort methylation-MAD filter, annotates retained sites with gene/TSS, GTF subfeature, cCRE, and CpG-island context, mean-imputes the limited remaining missing values per feature, produces raw beta-value and inverse-normal transformed phenotype BEDs, calculates phenotype PCs, and can merge those PCs with additional QTL covariates.

## What to provide

Use one pb-CpG-tools `.combined.bed.gz` file per sample. The workflow expects the `model` pileup output and uses `mod_score` as the methylation phenotype.

| pb-CpG-tools column | Workflow use |
| --- | --- |
| `#chrom`, `begin`, `end` | Site coordinates; retained as BED coordinates. |
| `mod_score` | Methylation value. It is multiplied by `0.01` by default to produce a 0–1 beta value. |
| `type` | Must have one value per input file. Use the sample's `.combined.bed.gz` output for standard meQTL analysis. |
| `cov` | Per-site coverage used for QC and site metadata. |
| `est_mod_count`, `est_unmod_count`, `discretized_mod_score` | Not loaded; these model-mode extras are not needed for QC or the QTL phenotype. |

All input BEDs must use the same reference genome and contig naming convention as the genotype data. The workflow removes contigs matching `X|Y|M|_` by default; set `FilterChroms` to an empty string to retain them.

## Sample manifest

Supply a TSV with a header and one row per sample:

```text
sample_id	file_path
SAMPLE_001	gs://my-bucket/SAMPLE_001.combined.bed.gz
SAMPLE_002	gs://my-bucket/SAMPLE_002.combined.bed.gz
```

The columns must be named `sample_id` and `file_path`. The workflow partitions the manifest into `SamplesPerShard` rows per task, and each task uses concurrent `gsutil` processes to localize only the `gs://` BEDs assigned to that shard. The Docker image includes compiled `crcmod`, allowing `gsutil` to use sliced downloads for large objects. Each shard writes one QC-flagged call table per autosome. Passing calls are derived from `per_sample_qc_pass` during the chromosome merge instead of uploading and localizing a duplicate filtered copy.

The workflow also requires a GTF, an ENCODE cCRE reference, and a UCSC `cpgIslandExt` reference built on the same genome assembly. The cCRE reference is a headerless six-column file: chromosome, BED start, BED end, V4 ID, V5 ID, and cCRE type. The cCRE intervals are interpreted as BED coordinates. This is deliberately broader than an enhancer annotation: V6 may also identify CTCF-only or DNase-H3K4me3 elements.

## Workflow design

```mermaid
flowchart LR
    A["Manifest with one BED per sample"] --> B["Split manifest into sample shards"]
    B --> C["Per-sample QC\nchromosome, coverage, extreme coverage"]
    C --> D["Autosome-split shard outputs"]
    D --> E["Parallel per-autosome cohort reductions\napply MinSampleFraction, MinSamples, and MAD"]
    E --> F["Header-aware streaming concatenation of metadata, QC, and BEDs"]
    F --> G["Gene/TSS, GTF subfeature, cCRE, and CpG-island\nannotation of passing sites"]
    G --> H["Phenotype PCs\noptional covariate merge"]
```

The cohort threshold is deliberately evaluated only after every parallel sample task has completed per-sample QC. This means the required number of samples is always:

```text
max(ceiling(total_samples × MinSampleFraction), MinSamples)
```

and never depends on parallel task boundaries. The threshold is evaluated separately within each autosome merge task, but every task uses the same complete-cohort denominator.

This chromosome-specific reduction also avoids the `data.table` limit of 2,147,483,647 rows that can be reached when all sample/site calls across the genome are combined into one long table. Each chromosome task localizes one call file per shard plus one small cohort-sample file; shard sample-QC files are localized only once by the final aggregation task.

## QC stages and run log

For every sample, the log reports:

1. Input site count and sites removed by the chromosome filter.
2. Sites failing `MinCoverage`.
3. Sites failing the extreme-coverage filter after meeting `MinCoverage`.
4. Sites passing both per-sample thresholds.

At cohort level, the log reports the union of sites after chromosome filtering, sites observed with adequate coverage in at least one sample, sites passing all per-sample QC in at least one sample, counts passing/failing the sample-presence threshold, counts failing the methylation-MAD threshold, and the number of sample/site values imputed.

The extreme-coverage threshold is a Tukey far-out fence calculated separately for each sample on `log10(cov)`. It is intended to exclude unusually high-coverage loci that may reflect copy-number or mapping artifacts.

## Important inputs

| Input | Default | Meaning |
| --- | --- | --- |
| `SamplesPerShard` | `25` | Number of samples processed by each parallel filtering task. |
| `ShardNumThreads` | `4` | Concurrent BED downloads per shard and threads available for shard parsing/decompression. |
| `MinCoverage` | `10` | Minimum `cov` required for a sample/site call. |
| `MinSampleFraction` | `0.95` | Fraction of the complete cohort that must pass per-sample QC for a site to be retained. Remaining QTL-BED missing values are imputed with the feature mean. |
| `MinSamples` | `0` | Optional additional minimum number of samples passing per-sample QC. |
| `MinMethylationMAD` | `0.003` | Minimum methylation MAD across per-sample-QC-passing observations required for QTL output. |
| `FenceK` | `3.0` | Far-out-fence multiplier used for extreme coverage. |
| `AutosomePrefix` | `chr` | Prefix used for autosome names in the input BEDs. Use `""` for BEDs with chromosomes named `1` through `22`. |
| `AnnotationGTF` | required | GTF gene model on the same reference build as the methylation calls. |
| `CCREAnnotations` | required | Headerless six-column ENCODE cCRE file: BED coordinates, V4 ID, V5 ID, and type. |
| `CpGIslandAnnotations` | required | UCSC `cpgIslandExt` table from the same reference build; accepts named columns or a headerless UCSC table. |
| `PromoterWindow` | `2000` | Bases either side of each strand-aware TSS classified as promoter. |
| `ValueColumn` | `mod_score` | Column used as methylation phenotype. |
| `ValueMultiplier` | `0.01` | Converts pb-CpG `mod_score` percentages to 0–1 beta values. |
| `AdditionalCovariates` | unset | Optional TSV containing `sample_id` plus genotype PCs or other covariates to merge with INT phenotype PCs. |
| `ShardMemoryGB` / `ShardDiskGB` | `64` / `250` | Resources for each parallel shard. Disk must accommodate the 25 input BEDs and chromosome-split outputs; `SamplesPerShard` remains 25 by default. |
| `MergeMemoryGB` / `MergeDiskGB` | `128` / `500` | Resources for each per-autosome cohort reduction and downstream phenotype-PC calculation. |
| `AggregateMemoryGB` / `AggregateDiskGB` | `64` / `1000` | Resources for final streaming aggregation. The larger disk accommodates simultaneous localization of all chromosome-level output families. |

## Outputs

| Output | Contents |
| --- | --- |
| `<prefix>.methylation.filtered.long.tsv.gz` | Calls from sites passing the final cohort threshold. |
| `<prefix>.methylation.site_qc.tsv.gz` | Compact all-site table with sample-presence counts and `keep_site`. |
| `<prefix>.methylation.site_metadata.tsv.gz` | All observed sites, including coverage and methylation means, standard deviations, CVs, methylation MAD, coverage fractions, sample counts, filter flags, `n_samples_imputed_in_qtl_bed`, and `keep_site`. |
| `<prefix>.methylation.sample_qc.tsv` | One row per sample with coverage filter counts, extreme-coverage cutoffs, and pass counts. |
| `<prefix>.methylation.filter_summary.tsv` | Counts of sites at each mutually exclusive cohort-QC stage. |
| `<prefix>.methylation.filter_counts.png` | Bar chart of the sequential cohort-QC counts. |
| `<prefix>.methylation.filter_upset.png` | UpSet-style chart showing overlap of low/missing coverage, extreme coverage, cohort sample-presence, and methylation-MAD conditions. |
| `<prefix>.methylation.passing_site_annotations.tsv.gz` | One row per retained site with nearest TSS/gene, promoter, GTF subfeatures, gene-body/intergenic, cCRE, and CpG-island annotations. |
| `<prefix>.methylation.raw.bed.gz` | TensorQTL-compatible phenotype BED with raw 0–1 methylation beta values. |
| `<prefix>.methylation.INT.bed.gz` | TensorQTL-compatible phenotype BED after site-wise rank-based inverse normal transformation. |
| `<prefix>.methylation_phenotype_PCs.INT.tsv` | Phenotype PCs calculated from the INT BED. |
| `<prefix>.methylation_QTL_covariates.INT.tsv` | Optional merged covariate matrix, written when `AdditionalCovariates` is supplied. |

The site metadata has two metric families:

- `*_all_calls`: all calls remaining after chromosome filtering, including calls that fail coverage QC.
- `*_passing_per_sample_qc`: only calls passing both the minimum- and extreme-coverage filters.

`fraction_samples_min_coverage` uses the complete input cohort as its denominator. `fraction_samples_passing_per_sample_qc` additionally excludes extreme-coverage calls. `pass_sample_presence_filter` and `pass_methylation_mad_filter` show the two cohort filters separately; `keep_site` is their final combined decision.

The filter-count chart and TSV use mutually exclusive stages so the counts add up to every observed site: insufficient minimum-coverage samples, loss of sufficient samples after extreme-coverage exclusions, low MAD after sample-presence QC, or passing all cohort filters. The UpSet-style plot is complementary: it retains overlapping conditions, including a site that still passes overall QC despite one or more missing or extreme-coverage calls.

## Gene, cCRE, and CpG-island annotation

The annotation task reads only `keep_site == TRUE` rows from the all-site metadata. It reports the nearest strand-aware gene TSS (`nearest_tss_distance`, gene ID/name, and strand), promoter overlaps within ±2 kb by default, gene-body overlaps, and a mutually exclusive `genomic_context` of `promoter`, `gene_body`, or `intergenic`. `gene_feature_context` further distinguishes `exon`, `intron`, and `gene_body_unclassified`; `in_exon`, `in_intron`, `in_cds`, `in_five_prime_utr`, `in_three_prime_utr`, and `in_utr` each have matching count and gene ID/name fields. Promoter takes priority when a site overlaps both a promoter and a gene body. Multiple overlapping genes are retained as semicolon-delimited IDs and names.

cCRE overlaps are reported separately because they may occur in any genic context. The annotation includes `in_ccre`, `n_overlapping_ccres`, `ccre_v4_id`, `ccre_v5_id`, and `ccre_type`; multiple overlaps are semicolon-delimited. The V4/V5 values are preserved exactly from the reference. `is_enhancer_like` identifies pELS/dELS overlaps, while `is_ctcf_only` identifies CTCF-only overlaps.

The UCSC `cpgIslandExt` table is read from its `chrom`, `chromStart`, `chromEnd`, and `name` columns (or their standard headerless column positions). `chromStart` and `chromEnd` are interpreted as UCSC BED coordinates. A passing CpG is classified with the following mutually exclusive priority:

1. `island`: the CpG overlaps an annotated CpG island.
2. `shore`: the CpG does not overlap an island but lies within 2,000 bp of an island boundary.
3. `shelf`: the CpG is not an island or shore site but lies 2,001–4,000 bp from an island boundary.
4. `open_sea`: the CpG is more than 4,000 bp from every annotated island, or its chromosome has no annotated island.

The output records this in `cpg_island_context`, together with `in_cpg_island`, overlapping island count/name, and the nearest island name and boundary distance. The classification is based on a 2 kb expansion of each island for shores and a 4 kb expansion for shelves; island overlap takes precedence, followed by shore, then shelf. All references must use the same genome build and chromosome naming convention as the PacBio calls.

## QTL phenotype BEDs, PCs, and covariates

Both phenotype BEDs use the first four columns required by TensorQTL:

```text
#chr	start	end	phenotype_id	SAMPLE_001	SAMPLE_002	...
chr1	10469	10470	chr1*10469*10470	0.73	0.68	...
```

The raw BED contains beta values (`mod_score / 100`). For every retained feature, samples missing after per-sample QC are filled with that feature's mean beta value among observed QC-passing samples. The metadata field `n_samples_imputed_in_qtl_bed` records exactly how many cells were imputed for that feature. The INT BED then rank-transforms each imputed CpG row across samples and is used for phenotype-PC calculation. The workflow calculates PCs only for the INT BED, following the existing molecular-QTL prepare workflows.

Set `AdditionalCovariates` to a TSV with a `sample_id` column to merge genotype PCs or other covariates with the INT phenotype PCs. The resulting covariate file has covariates as rows and sample IDs as columns, ready for TensorQTL.

The WDL defaults to `MinSampleFraction = 0.95` and `MinMethylationMAD = 0.003`. This retains sites present in at least 95% of samples and produces a complete matrix through feature-mean imputation before PCA and QTL mapping.

## Related files

- [`scripts/FilterMethylationShard.R`](../scripts/FilterMethylationShard.R): per-sample chromosome and coverage QC.
- [`scripts/MergeMethylationCohort.R`](../scripts/MergeMethylationCohort.R): cohort-wide filtering, metadata, imputation, phenotype BED, and QC plot generation.
- [`scripts/MethylationUtils.R`](../scripts/MethylationUtils.R): shared input, QC, transformation, and plotting functions sourced by both executable stages.
- [`scripts/AnnotateMethylationSites.R`](../scripts/AnnotateMethylationSites.R): passing-site gene/TSS, GTF-subfeature, cCRE, and CpG-island annotation task.
- [`workflows/merge_methylation.wdl`](../workflows/merge_methylation.wdl): WDL wrapper for parallel execution.
- [R script reference](scripts.md): reference for all project scripts.
- [Molecular QTL workflow reference](molecular-qtl-workflows.md): reference for all molecular workflow wrappers.
