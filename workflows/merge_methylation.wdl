version 1.0
import "calculate_phenotypePCs.wdl" as ComputePCs
import "MergeCovariates.wdl" as CovariateMerge

# The global merge is deliberately outside the scatter. A per-shard site
# filter would use a different denominator in every shard and would therefore
# not implement MinSampleFraction across the cohort.

task ShardMethylationManifest {
    input {
        File SampleManifest
        Int SamplesPerShard
    }

    command <<<
        set -euo pipefail

        # shellcheck disable=SC2016
        Rscript -e '
        library(data.table)
        manifest <- fread("~{SampleManifest}", sep = "\t", header = TRUE, quote = "", data.table = FALSE)
        required <- c("sample_id", "file_path")
        missing <- setdiff(required, names(manifest))
        if (length(missing) > 0) {
            stop("SampleManifest is missing required column(s): ", paste(missing, collapse = ", "))
        }
        manifest <- manifest[, required, drop = FALSE]
        if (nrow(manifest) < 1) stop("SampleManifest must contain at least one data row")
        if (anyNA(manifest$sample_id) || any(!nzchar(manifest$sample_id))) stop("SampleManifest contains an empty sample_id")
        if (any(!grepl("^[A-Za-z0-9._-]+$", manifest$sample_id))) {
            stop("sample_id values must match [A-Za-z0-9._-]+")
        }
        if (anyDuplicated(manifest$sample_id)) stop("Each sample_id must occur exactly once in SampleManifest")
        if (anyNA(manifest$file_path) || any(!nzchar(manifest$file_path))) stop("SampleManifest contains an empty file_path")
        shard_size <- as.integer("~{SamplesPerShard}")
        if (is.na(shard_size) || shard_size < 1) stop("SamplesPerShard must be at least 1")
        dir.create("shards", showWarnings = FALSE)
        starts <- seq.int(1L, nrow(manifest), by = shard_size)
        for (i in seq_along(starts)) {
            end <- min(nrow(manifest), starts[[i]] + shard_size - 1L)
            fwrite(manifest[starts[[i]]:end, , drop = FALSE],
                   sprintf("shards/methylation_manifest.shard.%05d.tsv", i - 1L), sep = "\t")
        }
        writeLines(as.character(nrow(manifest)), "total_samples.txt")
        '
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "2G"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        Array[File] ShardManifests = glob("shards/methylation_manifest.shard.*.tsv")
        Int TotalSamples = read_int("total_samples.txt")
    }
}

task FilterMethylationShard {
    input {
        File ShardManifest
        String OutputPrefix
        Float MinCoverage
        String FilterChroms
        Float FenceK
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail

        mkdir -p input_beds
        printf 'sample_id\tfile_path\n' > localized_manifest.tsv

        tail -n +2 "~{ShardManifest}" | while IFS=$'\t' read -r sample_id source_path; do
            [ -n "$sample_id" ] || continue
            local_path="input_beds/${sample_id}.combined.bed.gz"
            if [[ "$source_path" == gs://* ]]; then
                gsutil cp "$source_path" "$local_path"
            else
                if [ ! -f "$source_path" ]; then
                    echo "Input BED file for ${sample_id} is not accessible inside the task: ${source_path}" >&2
                    exit 1
                fi
                cp "$source_path" "$local_path"
            fi
            printf '%s\t%s\n' "$sample_id" "$local_path" >> localized_manifest.tsv
        done

        Rscript /tmp/FilterMethylationShard.R \
            --InputManifest localized_manifest.tsv \
            --OutputPrefix "~{OutputPrefix}" \
            --MinCoverage ~{MinCoverage} \
            --FilterChroms "~{FilterChroms}" \
            --FenceK ~{FenceK}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.per_sample_filtered.long.tsv.gz"
        File AllCalls = "~{OutputPrefix}.methylation.per_sample_qc.long.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
    }
}

task MergeMethylationShards {
    input {
        Array[File] FilteredCallShards
        Array[File] AllCallShards
        Array[File] SampleQcShards
        Int TotalSamples
        String OutputPrefix
        Float MinSampleFraction
        Int MinSamples
        Float MinMethylationMAD
        String ValueColumn
        Float ValueMultiplier
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' FilteredCallShards} > filtered_call_shards.list
        printf '%s\n' ~{sep=' ' AllCallShards} > all_call_shards.list
        printf '%s\n' ~{sep=' ' SampleQcShards} > sample_qc_shards.list

        Rscript /tmp/MergeMethylationCohort.R \
            --FilteredCallList filtered_call_shards.list \
            --AllCallList all_call_shards.list \
            --FilteredSampleQcList sample_qc_shards.list \
            --TotalSamples ~{TotalSamples} \
            --OutputPrefix "~{OutputPrefix}" \
            --MinSampleFraction ~{MinSampleFraction} \
            --MinSamples ~{MinSamples} \
            --MinMethylationMAD ~{MinMethylationMAD} \
            --ValueColumn "~{ValueColumn}" \
            --ValueMultiplier ~{ValueMultiplier}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.filtered.long.tsv.gz"
        File SiteQC = "~{OutputPrefix}.methylation.site_qc.tsv.gz"
        File SiteMetadata = "~{OutputPrefix}.methylation.site_metadata.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
        File FilterSummary = "~{OutputPrefix}.methylation.filter_summary.tsv"
        File FilterCountsPlot = "~{OutputPrefix}.methylation.filter_counts.png"
        File FilterUpsetPlot = "~{OutputPrefix}.methylation.filter_upset.png"
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
    }
}

task AnnotateMethylationSites {
    input {
        File SiteMetadata
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations
        String OutputPrefix
        Int PromoterWindow
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        Rscript /tmp/AnnotateMethylationSites.R \
            --SiteMetadata "~{SiteMetadata}" \
            --AnnotationGTF "~{AnnotationGTF}" \
            --CCREAnnotations "~{CCREAnnotations}" \
            --CpGIslandAnnotations "~{CpGIslandAnnotations}" \
            --OutputPrefix "~{OutputPrefix}" \
            --PromoterWindow ~{PromoterWindow}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File PassingSiteAnnotations = "~{OutputPrefix}.methylation.passing_site_annotations.tsv.gz"
    }
}

workflow MergeMethylation {
    input {
        # TSV with sample_id and file_path columns. gs:// files are localized
        # inside each shard task with gsutil.
        File SampleManifest
        String OutputPrefix
        File? AdditionalCovariates
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations

        Int SamplesPerShard = 25
        Float MinCoverage = 10.0
        Float MinSampleFraction = 0.95
        Int MinSamples = 0
        Float MinMethylationMAD = 0.003
        String FilterChroms = "X|Y|M|_"
        Float FenceK = 3.0
        Int PromoterWindow = 2000
        String ValueColumn = "mod_score"
        Float ValueMultiplier = 0.01

        Int ShardMemoryGB = 16
        Int ShardDiskGB = 100
        Int MergeMemoryGB = 64
        Int MergeDiskGB = 200
        Int AnnotationMemoryGB = 16
        Int AnnotationDiskGB = 100
        Int NumThreads = 1
    }

    call ShardMethylationManifest {
        input:
            SampleManifest = SampleManifest,
            SamplesPerShard = SamplesPerShard
    }

    scatter (shard_index in range(length(ShardMethylationManifest.ShardManifests))) {
        File shard_manifest = ShardMethylationManifest.ShardManifests[shard_index]
        String shard_output_prefix = "~{OutputPrefix}.shard.~{shard_index}"

        call FilterMethylationShard {
            input:
                ShardManifest = shard_manifest,
                OutputPrefix = shard_output_prefix,
                MinCoverage = MinCoverage,
                FilterChroms = FilterChroms,
                FenceK = FenceK,
                MemoryGB = ShardMemoryGB,
                DiskGB = ShardDiskGB,
                NumThreads = NumThreads
        }
    }

    call MergeMethylationShards {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCalls,
            AllCallShards = FilterMethylationShard.AllCalls,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            OutputPrefix = OutputPrefix,
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call AnnotateMethylationSites {
        input:
            SiteMetadata = MergeMethylationShards.SiteMetadata,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            OutputPrefix = OutputPrefix,
            PromoterWindow = PromoterWindow,
            MemoryGB = AnnotationMemoryGB,
            DiskGB = AnnotationDiskGB
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = MergeMethylationShards.IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = MergeMemoryGB,
            disk_space = MergeDiskGB,
            num_threads = NumThreads
    }

    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT"
        }
    }

    output {
        File FilteredCalls = MergeMethylationShards.FilteredCalls
        File SiteQC = MergeMethylationShards.SiteQC
        File SiteMetadata = MergeMethylationShards.SiteMetadata
        File SampleQC = MergeMethylationShards.SampleQC
        File FilterSummary = MergeMethylationShards.FilterSummary
        File FilterCountsPlot = MergeMethylationShards.FilterCountsPlot
        File FilterUpsetPlot = MergeMethylationShards.FilterUpsetPlot
        File RawMethylationBed = MergeMethylationShards.RawMethylationBed
        File IntMethylationBed = MergeMethylationShards.IntMethylationBed
        File PassingSiteAnnotations = AnnotateMethylationSites.PassingSiteAnnotations
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        Array[File] ShardSampleQC = FilterMethylationShard.SampleQC
    }
}
