version 1.0

struct PlinkChromosome {
    String chromosome
    String format
    File genotype
    File variants
    File samples
}

workflow TransLDRegions {
    input {
        Array[PlinkChromosome] genotypes
        Array[File] association_files
        File ancestry_samples
        File chrom_sizes
        String genome_build
        Float pvalue_threshold
        String trans_ld_prepare_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        String trans_ld_ancestry_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        String trans_ld_combine_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        String variant_column = "variant_id"
        String phenotype_column = "phenotype_id"
        String pvalue_column = "pval"
        String chrom_column = "chrom"
        String position_column = "pos"
        String modality_column = "modality"
        String dataset_column = "dataset_id"
        String cell_type_column = "cell_type"
        String default_modality = "unspecified"
        Float maf_threshold = 0.01
        Float r2_threshold = 0.1
        Int initial_search_bp = 1000000
        Int max_search_bp = 10000000
        Int edge_bp = 100000
        Int padding_bp = 100000
        Int fallback_bp = 1000000
        Int min_region_bp = 2000000
        Int threads = 4
        Int ld_memory_gb = 16
        Int ld_disk_gb = 100
        Int preparation_memory_gb = 8
        Int preparation_disk_gb = 50
        Int merge_memory_gb = 16
        Int merge_disk_gb = 50
        Int preemptible_attempts = 1
    }
    scatter (shard in genotypes) {
        String chromosome_name = shard.chromosome
    }
    call PrepareInputs {
        input:
            association_files = association_files,
            ancestry_samples = ancestry_samples,
            chrom_sizes = chrom_sizes,
            chromosome_names = chromosome_name,
            genome_build = genome_build,
            pvalue_threshold = pvalue_threshold,
            variant_column = variant_column,
            phenotype_column = phenotype_column,
            pvalue_column = pvalue_column,
            chrom_column = chrom_column,
            position_column = position_column,
            modality_column = modality_column,
            dataset_column = dataset_column,
            cell_type_column = cell_type_column,
            default_modality = default_modality,
            docker_image = trans_ld_prepare_image,
            memory_gb = preparation_memory_gb,
            disk_gb = preparation_disk_gb,
            preemptible_attempts = preemptible_attempts
    }
    scatter (shard in genotypes) {
        scatter (group in read_lines(PrepareInputs.ancestries)) {
            call AncestryLD {
                input:
                    genotype = shard.genotype,
                    variants = shard.variants,
                    samples = shard.samples,
                    genotype_format = shard.format,
                    chromosome = shard.chromosome,
                    ancestry_samples = ancestry_samples,
                    ancestry_group = group,
                    associations = PrepareInputs.associations,
                    chrom_sizes = chrom_sizes,
                    maf_threshold = maf_threshold,
                    r2_threshold = r2_threshold,
                    initial_search_bp = initial_search_bp,
                    max_search_bp = max_search_bp,
                    edge_bp = edge_bp,
                    threads = threads,
                    memory_gb = ld_memory_gb,
                    disk_gb = ld_disk_gb,
                    docker_image = trans_ld_ancestry_image,
                    preemptible_attempts = preemptible_attempts
            }
        }
    }
    call CombineRegions {
        input:
            associations = PrepareInputs.associations,
            interval_files = flatten(AncestryLD.intervals),
            chrom_sizes = chrom_sizes,
            genome_build = genome_build,
            padding_bp = padding_bp,
            fallback_bp = fallback_bp,
            min_region_bp = min_region_bp,
            docker_image = trans_ld_combine_image,
            memory_gb = merge_memory_gb,
            disk_gb = merge_disk_gb,
            preemptible_attempts = preemptible_attempts
    }
    output {
        File regions = CombineRegions.regions
        File regions_bed = CombineRegions.regions_bed
        File region_associations = CombineRegions.region_associations
        File trans_window_associations = CombineRegions.trans_window_associations
        File seed_regions = CombineRegions.seed_regions
        File ancestry_intervals = CombineRegions.ancestry_intervals
        File region_qc = CombineRegions.qc
        File preparation_qc = PrepareInputs.qc
        File normalized_associations = PrepareInputs.associations
        File ancestry_groups = PrepareInputs.ancestries
        Array[File] ancestry_qc = flatten(AncestryLD.qc)
        Array[File] ancestry_allele_counts = flatten(AncestryLD.allele_counts)
        Array[File] ld_logs = flatten(AncestryLD.log_file)
        File preparation_log = PrepareInputs.log_file
        File merge_log = CombineRegions.log_file
    }
    parameter_meta {
        genotypes: "One entry per chromosome (1-22 or X). format is bed or pgen. All three genotype components must be explicit File values."
        association_files: "TSV, TSV.gz or Parquet. Use selected trans associations; pvalue_threshold only filters rows and does not calculate significance."
        ancestry_samples: "TSV with header IID and ancestry; unique IID, one group per participant."
        chrom_sizes: "Two-column chromosome/length file or FASTA .fai, for the specified genome build."
        trans_ld_prepare_image: "Pinned image for input preparation."
        trans_ld_ancestry_image: "Pinned image for ancestry LD."
        trans_ld_combine_image: "Pinned image for region combination."
        maf_threshold: "Strict lower bound on within-group minor allele frequency. No MAC cutoff."
        min_region_bp: "Minimum final interval width before overlap merging; default 2 Mb total."
        max_search_bp: "Maximum LD search distance per side. Search-limit hits are retained and flagged."
    }
}

task PrepareInputs {
    input {
        Array[File] association_files
        File ancestry_samples
        File chrom_sizes
        Array[String] chromosome_names
        String genome_build
        Float pvalue_threshold
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        String variant_column = "variant_id"
        String phenotype_column = "phenotype_id"
        String pvalue_column = "pval"
        String chrom_column = "chrom"
        String position_column = "pos"
        String modality_column = "modality"
        String dataset_column = "dataset_id"
        String cell_type_column = "cell_type"
        String default_modality = "unspecified"
        Int memory_gb = 8
        Int disk_gb = 50
        Int preemptible_attempts = 1
    }
command <<<
        set -euo pipefail
        echo '[prepare] Starting association and ancestry validation'
        python3 -m trans_ld_regions.prepare \
            --association-list '~{write_lines(association_files)}' \
            --ancestry '~{sub(ancestry_samples, "'", "'\"'\"'")}' \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --chromosomes '~{write_lines(chromosome_names)}' \
            --genome-build '~{sub(genome_build, "'", "'\"'\"'")}' \
            --pvalue-threshold ~{pvalue_threshold} \
            --variant-column '~{sub(variant_column, "'", "'\"'\"'")}' \
            --phenotype-column '~{sub(phenotype_column, "'", "'\"'\"'")}' \
            --pvalue-column '~{sub(pvalue_column, "'", "'\"'\"'")}' \
            --chrom-column '~{sub(chrom_column, "'", "'\"'\"'")}' \
            --position-column '~{sub(position_column, "'", "'\"'\"'")}' \
            --modality-column '~{sub(modality_column, "'", "'\"'\"'")}' \
            --dataset-column '~{sub(dataset_column, "'", "'\"'\"'")}' \
            --cell-type-column '~{sub(cell_type_column, "'", "'\"'\"'")}' \
            --default-modality '~{sub(default_modality, "'", "'\"'\"'")}' \
            2>&1 | tee preparation.log
        echo '[prepare] Complete'
    >>>
    output {
        File associations = "associations.tsv.gz"
        File ancestries = "ancestries.txt"
        File qc = "preparation.json"
        File log_file = "preparation.log"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible_attempts
    }
}

task AncestryLD {
    input {
        File genotype
        File variants
        File samples
        String genotype_format
        String chromosome
        File ancestry_samples
        String ancestry_group
        File associations
        File chrom_sizes
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        Float maf_threshold = 0.01
        Float r2_threshold = 0.1
        Int initial_search_bp = 1000000
        Int max_search_bp = 10000000
        Int edge_bp = 100000
        Int threads = 4
        Int memory_gb = 16
        Int disk_gb = 100
        Int preemptible_attempts = 1
    }
command <<<
        set -euo pipefail
        echo '[ld] Starting ancestry-specific LD analysis'
        python3 -m trans_ld_regions.ld \
            --genotype '~{sub(genotype, "'", "'\"'\"'")}' \
            --variants '~{sub(variants, "'", "'\"'\"'")}' \
            --samples '~{sub(samples, "'", "'\"'\"'")}' \
            --format '~{sub(genotype_format, "'", "'\"'\"'")}' \
            --chromosome '~{sub(chromosome, "'", "'\"'\"'")}' \
            --ancestry '~{sub(ancestry_samples, "'", "'\"'\"'")}' \
            --group '~{sub(ancestry_group, "'", "'\"'\"'")}' \
            --associations '~{sub(associations, "'", "'\"'\"'")}' \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --maf-threshold ~{maf_threshold} --r2-threshold ~{r2_threshold} \
            --initial-search-bp ~{initial_search_bp} --max-search-bp ~{max_search_bp} \
            --edge-bp ~{edge_bp} --threads ~{threads} --memory-mb ~{memory_gb * 800} \
            2>&1 | tee ancestry_ld.log
        echo '[ld] Complete'
    >>>
    output {
        File intervals = "intervals.tsv"
        File allele_counts = "allele_counts.tsv.gz"
        File qc = "ld_qc.json"
        File log_file = "ancestry_ld.log"
    }
    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible_attempts
    }
}

task CombineRegions {
    input {
        File associations
        Array[File] interval_files
        File chrom_sizes
        String genome_build
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-trans-ld-regions@sha256:523c975d2b616b2b2f894e9ee04313174708c8dce84c82cd671c1f8a2064d280"
        Int padding_bp = 100000
        Int fallback_bp = 1000000
        Int min_region_bp = 2000000
        Int memory_gb = 16
        Int disk_gb = 50
        Int preemptible_attempts = 1
    }
command <<<
        set -euo pipefail
        echo '[merge] Combining ancestry boundaries and enforcing minimum width'
        python3 -m trans_ld_regions.combine \
            --associations '~{sub(associations, "'", "'\"'\"'")}' \
            --interval-list '~{write_lines(interval_files)}' \
            --chrom-sizes '~{sub(chrom_sizes, "'", "'\"'\"'")}' \
            --genome-build '~{sub(genome_build, "'", "'\"'\"'")}' \
            --padding-bp ~{padding_bp} --fallback-bp ~{fallback_bp} \
            --min-region-bp ~{min_region_bp} \
            2>&1 | tee merge.log
        echo '[merge] Complete'
    >>>
    output {
        File regions = "regions.tsv"
        File regions_bed = "regions.bed"
        File region_associations = "region_associations.tsv.gz"
        File trans_window_associations = "trans_window_associations.tsv.gz"
        File seed_regions = "seed_regions.tsv"
        File ancestry_intervals = "ancestry_intervals.tsv.gz"
        File qc = "region_qc.json"
        File log_file = "merge.log"
    }
    runtime {
        docker: docker_image
        cpu: 1
        memory: "~{memory_gb} GiB"
        disks: "local-disk ~{disk_gb} HDD"
        preemptible: preemptible_attempts
    }
}
