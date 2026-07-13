# PacBio 5mC QTL workflow

[Back to main README](../README.md)

This guide describes how to prepare pb-CpG-tools site-level 5mC calls for molecular QTL mapping. The workflow applies coverage QC per sample, filters sites against the whole cohort, applies a cohort methylation-MAD filter, reports a diagnostic methylation–coverage correlation, annotates retained sites with gene/TSS, GTF subfeature, cCRE, and CpG-island context, mean-imputes the limited remaining missing values per feature, produces raw beta-value and inverse-normal transformed phenotype BEDs, calculates phenotype PCs, and can merge those PCs with additional QTL covariates.

## What to provide

Use one pb-CpG-tools `.combined.bed.gz` file per sample. The workflow expects the `model` pileup output and uses `mod_score` as the methylation phenotype.

## Workflow entry points

Choose the entry point that matches how inputs are organized:

- [`workflows/methylation/ProcessMethylationSample.wdl`](../workflows/methylation/ProcessMethylationSample.wdl) is the Terra-table sample workflow. Run it once per sample entity with `SampleID` and `MethylationBed` bound directly to table columns. It has no manifest input and writes one sample-QC file plus one QC-flagged call table for each autosome. Its per-sample parsing/filtering is performed by the compiled streaming Rust tool in `rust/methylation_filter`, using the slim `prepare_qtl-methylation-rust` image. The task reserves four CPUs: one parses the input and three perform bounded parallel chromosome-output compression. The manifest-sharded workflow continues to use the established R implementation.
- [`workflows/methylation/AggregateMethylationCohort.wdl`](../workflows/methylation/AggregateMethylationCohort.wdl) is the Terra cohort workflow. Supply one `CohortManifest` file containing the per-sample QC and 22 chromosome-output paths. This avoids placing 22 large file arrays into Terra's workflow-submission JSON. Its per-chromosome merge uses a bounded-memory Rust k-way merge, including staged intermediate runs when more than 128 sample files must be opened.
- [`workflows/methylation/merge_methylation.wdl`](../workflows/methylation/merge_methylation.wdl) remains the manifest/shard entry point. It retains its array-based internal aggregation implementation, which is safe because those file arrays are created inside Cromwell rather than submitted through Terra's API.

The cohort manifest must be a tab-separated file with exactly these columns. Each path should be a durable `gs://` path from the corresponding `ProcessMethylationSample` output:

```text
sample_id	sample_qc	autosome01	autosome02	...	autosome22
1000234	gs://.../1000234.methylation.sample_qc.tsv	gs://.../1000234.methylation.autosome01.per_sample_qc.long.tsv.gz	gs://.../1000234.methylation.autosome02.per_sample_qc.long.tsv.gz	...	gs://.../1000234.methylation.autosome22.per_sample_qc.long.tsv.gz
```

The workflow validates and splits this manifest internally, localizing only the chromosome call files needed by each chromosome merge. Sample-QC files are consolidated once before the chromosome scatter.

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

This chromosome-specific reduction also avoids the `data.table` limit of 2,147,483,647 rows that can be reached when all sample/site calls across the genome are combined into one long table. Each chromosome task localizes one call file per shard plus one consolidated cohort sample-QC table; the original shard sample-QC files are localized only once to build that consolidated table.

## Local CpG-correlation clustering and sample connectivity

[`scripts/methylation/AnalyzeMethylationCpGCorrelation.R`](../scripts/methylation/AnalyzeMethylationCpGCorrelation.R) runs independently on each chromosome INT BED after preliminary phenotype PCs and optional additional covariates have been assembled. It residualizes CpGs against that full covariate set, calls local correlation clusters, and selects the CpG with the greatest sum of absolute local correlations in each cluster. Unclustered CpGs are retained as singleton representatives.

```bash
Rscript AnalyzeMethylationCpGCorrelation.R \
  --InputBed cohort.methylation.INT.bed.gz \
  --Covariates cohort.methylation_QTL_covariates.INT.tsv \
  --OutputPrefix cohort \
  --WindowBP 1000 \
  --MinAbsCorrelation 0.95
```

The default is a conservative `|r| >= 0.95` threshold. The workflow writes per-chromosome cluster tables, summaries, and QC plots, then uses every representative CpG for landmark-based sample connectivity by default. These representatives are used only for connectivity; correlated CpGs are not removed from either final BED. Samples with connectivity Z-score below −3 are removed consistently from the final raw BED, INT BED, filtered long calls, and sample-QC output. Final phenotype PCs and optional QTL covariates are recalculated after this sample filter.

## QC stages and run log

For every sample, the log reports:

1. Input site count and sites removed by the chromosome filter.
2. Sites failing `MinCoverage`.
3. Sites failing the extreme-coverage filter after meeting `MinCoverage`.
4. Sites passing both per-sample thresholds.

At cohort level, the log reports the union of sites after chromosome filtering, sites observed with adequate coverage in at least one sample, sites passing all per-sample QC in at least one sample, counts passing/failing the sample-presence threshold, counts failing the methylation-MAD threshold, and the number of sample/site values imputed.

The extreme-coverage threshold is a Tukey far-out fence calculated separately for each sample on `log10(cov)`. It is intended to exclude unusually high-coverage loci that may reflect copy-number or mapping artifacts.

During each chromosome merge, the cohort-metric reduction reports progress every 1,000 CpGs by default. The standalone `MergeMethylationCohort.R` option `--ProgressEverySites` can adjust that interval.

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
| `CorrelationWindowBP` / `CorrelationMinAbsCorrelation` | `1000` / `0.95` | Local window and absolute Pearson threshold used to form residualized CpG clusters. |
| `MaxConnectivityFeatures` / `ConnectivityLandmarks` / `ConnectivityZThreshold` | `0` / `200` / `-3` | Optional representative-CpG cap (`0` uses all representatives), landmark sample count, and low-connectivity outlier threshold. |

## Outputs

| Output | Contents |
| --- | --- |
| `<prefix>.methylation.filtered.long.tsv.gz` | Calls from sites passing cohort QC and the sample connectivity filter. |
| `<prefix>.methylation.site_qc.tsv.gz` | Compact all-site table with sample-presence counts and `keep_site`. |
| `<prefix>.methylation.site_metadata.tsv.gz` | All observed sites, including coverage and methylation means, standard deviations, CVs, methylation MAD, coverage fractions, sample counts, `coverage_methylation_spearman_rho`, filter flags, `n_samples_imputed_in_qtl_bed`, and `keep_site`. |
| `<prefix>.methylation.passing_site_metadata.tsv.gz` | Streamed subset of the site metadata containing only `keep_site == TRUE`; consumed by the annotation task. |
| `<prefix>.methylation.sample_qc.tsv` | One row per sample with coverage QC plus connectivity score, Z-score, and pass flag. |
| `<prefix>.methylation.filter_summary.tsv` | Counts of sites at each mutually exclusive cohort-QC stage. |
| `<prefix>.methylation.filter_counts.png` | Bar chart of the sequential cohort-QC counts. |
| `<prefix>.methylation.filter_upset.png` | ggupset UpSet chart showing overlap of low/missing coverage, extreme coverage, cohort sample-presence, and methylation-MAD conditions. |
| `<prefix>.methylation.passing_site_annotations.tsv.gz` | One row per retained site with nearest TSS/gene, promoter, GTF subfeatures, gene-body/intergenic, cCRE, and CpG-island annotations. |
| `<prefix>.methylation.raw.bed.gz` | TensorQTL-compatible raw beta-value BED after sample connectivity filtering. |
| `<prefix>.methylation.INT.bed.gz` | TensorQTL-compatible INT BED after the same sample connectivity filtering. |
| `<prefix>.methylation.connectivity_representative_cpgs.tsv.gz` | Selected cluster representatives used for connectivity, including local-connectivity values. |
| `<prefix>.methylation.connectivity_outliers.tsv` / `.connectivity_summary.tsv` | Removed samples and connectivity method/threshold summary. |
| `<prefix>.methylation_phenotype_PCs.INT.tsv` | Phenotype PCs calculated from the INT BED. |
| `<prefix>.methylation_QTL_covariates.INT.tsv` | Optional merged covariate matrix, written when `AdditionalCovariates` is supplied. |

The site metadata has two metric families:

- `*_all_calls`: all calls remaining after chromosome filtering, including calls that fail coverage QC.
- `*_passing_per_sample_qc`: only calls passing both the minimum- and extreme-coverage filters.

`fraction_samples_min_coverage` uses the complete input cohort as its denominator. `fraction_samples_passing_per_sample_qc` additionally excludes extreme-coverage calls. `pass_sample_presence_filter` and `pass_methylation_mad_filter` show the two cohort filters separately; `keep_site` is their final combined decision.

`coverage_methylation_spearman_rho` is a diagnostic only: it is the per-site Spearman correlation between methylation and `log1p(cov) - log1p(sample_median_cov)`, using calls that passed the per-sample coverage QC. `n_samples_coverage_methylation_correlation` gives the number of calls contributing to that statistic. Neither field changes site retention or QTL output.

The filter-count chart and TSV use mutually exclusive stages so the counts add up to every observed site: insufficient minimum-coverage samples, loss of sufficient samples after extreme-coverage exclusions, low MAD after sample-presence QC, or passing all cohort filters. The ggupset UpSet plot is complementary: it retains overlapping conditions, including a site that still passes overall QC despite one or more missing or extreme-coverage calls.

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

Set `AdditionalCovariates` to a TSV with a `sample_id` column to merge genotype PCs or other covariates with the INT phenotype PCs. The file may contain samples outside the methylation cohort; the workflow discards those rows automatically. Every methylation sample must have a matching row in `AdditionalCovariates`. The resulting covariate file has covariates as rows and sample IDs as columns, ready for TensorQTL.

The WDL defaults to `MinSampleFraction = 0.95` and `MinMethylationMAD = 0.003`. This retains sites present in at least 95% of samples and produces a complete matrix through feature-mean imputation before PCA and QTL mapping.

## Related files

- [`scripts/methylation/FilterMethylationShard.R`](../scripts/methylation/FilterMethylationShard.R): per-sample chromosome and coverage QC.
- [`scripts/methylation/BuildMethylationCohortSamples.R`](../scripts/methylation/BuildMethylationCohortSamples.R): builds the cohort sample order and one consolidated sample-QC table from sample or shard QC outputs.
- [`scripts/methylation/AnalyzeMethylationCpGCorrelation.R`](../scripts/methylation/AnalyzeMethylationCpGCorrelation.R): calls covariate-adjusted local CpG clusters and selects the most connected representative per cluster.
- [`scripts/methylation/BuildMethylationCorrelationCovariates.R`](../scripts/methylation/BuildMethylationCorrelationCovariates.R): formats preliminary phenotype PCs plus optional additional covariates for correlation clustering.
- [`scripts/methylation/FinalizeMethylationConnectivity.R`](../scripts/methylation/FinalizeMethylationConnectivity.R): calculates landmark sample connectivity from representatives and filters raw/INT BEDs plus long calls to passing samples.
- [`scripts/methylation/MergeMethylationCohort.R`](../scripts/methylation/MergeMethylationCohort.R): cohort-wide filtering, metadata, imputation, phenotype BED, and QC plot generation.
- [`scripts/methylation/MethylationUtils.R`](../scripts/methylation/MethylationUtils.R): shared input, QC, transformation, and plotting functions sourced by both executable stages.
- [`scripts/methylation/AnnotateMethylationSites.R`](../scripts/methylation/AnnotateMethylationSites.R): passing-site gene/TSS, GTF-subfeature, cCRE, and CpG-island annotation task.
- [`workflows/methylation/ProcessMethylationSample.wdl`](../workflows/methylation/ProcessMethylationSample.wdl), [`workflows/methylation/AggregateMethylationCohort.wdl`](../workflows/methylation/AggregateMethylationCohort.wdl), and [`workflows/methylation/merge_methylation.wdl`](../workflows/methylation/merge_methylation.wdl): sample-table, cohort, and manifest/shard WDL entry points.
- [`workflows/methylation/cohort_aggregation.wdl`](../workflows/methylation/cohort_aggregation.wdl), [`workflows/methylation/connectivity.wdl`](../workflows/methylation/connectivity.wdl), [`workflows/methylation/annotation.wdl`](../workflows/methylation/annotation.wdl), and [`workflows/methylation/qtl_covariates.wdl`](../workflows/methylation/qtl_covariates.wdl): internal stage workflows called by the cohort entry point for aggregation, connectivity analysis, site annotation, and final QTL covariate preparation, respectively.
- [R script reference](scripts.md): reference for all project scripts.
- [Molecular QTL workflow reference](molecular-qtl-workflows.md): reference for all molecular workflow wrappers.
