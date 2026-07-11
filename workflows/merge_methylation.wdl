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
        String AutosomePrefix
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
            --FenceK ~{FenceK} \
            --AutosomePrefix "~{AutosomePrefix}"
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
        File FilteredCallsAutosome01 = "~{OutputPrefix}.methylation.autosome01.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome02 = "~{OutputPrefix}.methylation.autosome02.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome03 = "~{OutputPrefix}.methylation.autosome03.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome04 = "~{OutputPrefix}.methylation.autosome04.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome05 = "~{OutputPrefix}.methylation.autosome05.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome06 = "~{OutputPrefix}.methylation.autosome06.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome07 = "~{OutputPrefix}.methylation.autosome07.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome08 = "~{OutputPrefix}.methylation.autosome08.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome09 = "~{OutputPrefix}.methylation.autosome09.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome10 = "~{OutputPrefix}.methylation.autosome10.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome11 = "~{OutputPrefix}.methylation.autosome11.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome12 = "~{OutputPrefix}.methylation.autosome12.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome13 = "~{OutputPrefix}.methylation.autosome13.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome14 = "~{OutputPrefix}.methylation.autosome14.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome15 = "~{OutputPrefix}.methylation.autosome15.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome16 = "~{OutputPrefix}.methylation.autosome16.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome17 = "~{OutputPrefix}.methylation.autosome17.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome18 = "~{OutputPrefix}.methylation.autosome18.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome19 = "~{OutputPrefix}.methylation.autosome19.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome20 = "~{OutputPrefix}.methylation.autosome20.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome21 = "~{OutputPrefix}.methylation.autosome21.per_sample_filtered.long.tsv.gz"
        File FilteredCallsAutosome22 = "~{OutputPrefix}.methylation.autosome22.per_sample_filtered.long.tsv.gz"
        File AllCallsAutosome01 = "~{OutputPrefix}.methylation.autosome01.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome02 = "~{OutputPrefix}.methylation.autosome02.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome03 = "~{OutputPrefix}.methylation.autosome03.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome04 = "~{OutputPrefix}.methylation.autosome04.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome05 = "~{OutputPrefix}.methylation.autosome05.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome06 = "~{OutputPrefix}.methylation.autosome06.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome07 = "~{OutputPrefix}.methylation.autosome07.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome08 = "~{OutputPrefix}.methylation.autosome08.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome09 = "~{OutputPrefix}.methylation.autosome09.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome10 = "~{OutputPrefix}.methylation.autosome10.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome11 = "~{OutputPrefix}.methylation.autosome11.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome12 = "~{OutputPrefix}.methylation.autosome12.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome13 = "~{OutputPrefix}.methylation.autosome13.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome14 = "~{OutputPrefix}.methylation.autosome14.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome15 = "~{OutputPrefix}.methylation.autosome15.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome16 = "~{OutputPrefix}.methylation.autosome16.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome17 = "~{OutputPrefix}.methylation.autosome17.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome18 = "~{OutputPrefix}.methylation.autosome18.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome19 = "~{OutputPrefix}.methylation.autosome19.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome20 = "~{OutputPrefix}.methylation.autosome20.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome21 = "~{OutputPrefix}.methylation.autosome21.per_sample_qc.long.tsv.gz"
        File AllCallsAutosome22 = "~{OutputPrefix}.methylation.autosome22.per_sample_qc.long.tsv.gz"
    }
}

task MergeMethylationChromosome {
    input {
        Array[File] FilteredCallShards
        Array[File] AllCallShards
        Array[File] SampleQcShards
        Int TotalSamples
        String Chromosome
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
            --Chromosome "~{Chromosome}" \
            --OutputPrefix "~{OutputPrefix}" \
            --MinSampleFraction ~{MinSampleFraction} \
            --MinSamples ~{MinSamples} \
            --MinMethylationMAD ~{MinMethylationMAD} \
            --ValueColumn "~{ValueColumn}" \
            --ValueMultiplier ~{ValueMultiplier} \
            --SkipSampleQC \
            --SkipFilterPlots
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
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
    }
}

task AggregateMethylationChromosomes {
    input {
        Array[File] FilteredCallsByChromosome
        Array[File] SiteQCByChromosome
        Array[File] SiteMetadataByChromosome
        Array[File] RawMethylationBedByChromosome
        Array[File] IntMethylationBedByChromosome
        Array[File] SampleQcShards
        Int TotalSamples
        String OutputPrefix
        Int MemoryGB
        Int DiskGB
        Int NumThreads
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' FilteredCallsByChromosome} > filtered_calls_by_chromosome.list
        printf '%s\n' ~{sep=' ' SiteQCByChromosome} > site_qc_by_chromosome.list
        printf '%s\n' ~{sep=' ' SiteMetadataByChromosome} > site_metadata_by_chromosome.list
        printf '%s\n' ~{sep=' ' RawMethylationBedByChromosome} > raw_beds_by_chromosome.list
        printf '%s\n' ~{sep=' ' IntMethylationBedByChromosome} > int_beds_by_chromosome.list
        printf '%s\n' ~{sep=' ' SampleQcShards} > sample_qc_shards.list

        Rscript /tmp/AggregateMethylationChromosomes.R \
            --FilteredCallList filtered_calls_by_chromosome.list \
            --SiteQcList site_qc_by_chromosome.list \
            --SiteMetadataList site_metadata_by_chromosome.list \
            --RawBedList raw_beds_by_chromosome.list \
            --IntBedList int_beds_by_chromosome.list \
            --SampleQcList sample_qc_shards.list \
            --TotalSamples ~{TotalSamples} \
            --OutputPrefix "~{OutputPrefix}"
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
        String AutosomePrefix = "chr"
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
                AutosomePrefix = AutosomePrefix,
                MemoryGB = ShardMemoryGB,
                DiskGB = ShardDiskGB,
                NumThreads = NumThreads
        }
    }

    call MergeMethylationChromosome as MergeMethylationAutosome01 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome01,
            AllCallShards = FilterMethylationShard.AllCallsAutosome01,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "1",
            OutputPrefix = OutputPrefix + ".autosome01",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome02 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome02,
            AllCallShards = FilterMethylationShard.AllCallsAutosome02,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "2",
            OutputPrefix = OutputPrefix + ".autosome02",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome03 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome03,
            AllCallShards = FilterMethylationShard.AllCallsAutosome03,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "3",
            OutputPrefix = OutputPrefix + ".autosome03",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome04 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome04,
            AllCallShards = FilterMethylationShard.AllCallsAutosome04,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "4",
            OutputPrefix = OutputPrefix + ".autosome04",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome05 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome05,
            AllCallShards = FilterMethylationShard.AllCallsAutosome05,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "5",
            OutputPrefix = OutputPrefix + ".autosome05",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome06 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome06,
            AllCallShards = FilterMethylationShard.AllCallsAutosome06,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "6",
            OutputPrefix = OutputPrefix + ".autosome06",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome07 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome07,
            AllCallShards = FilterMethylationShard.AllCallsAutosome07,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "7",
            OutputPrefix = OutputPrefix + ".autosome07",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome08 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome08,
            AllCallShards = FilterMethylationShard.AllCallsAutosome08,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "8",
            OutputPrefix = OutputPrefix + ".autosome08",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome09 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome09,
            AllCallShards = FilterMethylationShard.AllCallsAutosome09,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "9",
            OutputPrefix = OutputPrefix + ".autosome09",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome10 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome10,
            AllCallShards = FilterMethylationShard.AllCallsAutosome10,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "10",
            OutputPrefix = OutputPrefix + ".autosome10",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome11 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome11,
            AllCallShards = FilterMethylationShard.AllCallsAutosome11,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "11",
            OutputPrefix = OutputPrefix + ".autosome11",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome12 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome12,
            AllCallShards = FilterMethylationShard.AllCallsAutosome12,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "12",
            OutputPrefix = OutputPrefix + ".autosome12",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome13 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome13,
            AllCallShards = FilterMethylationShard.AllCallsAutosome13,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "13",
            OutputPrefix = OutputPrefix + ".autosome13",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome14 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome14,
            AllCallShards = FilterMethylationShard.AllCallsAutosome14,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "14",
            OutputPrefix = OutputPrefix + ".autosome14",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome15 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome15,
            AllCallShards = FilterMethylationShard.AllCallsAutosome15,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "15",
            OutputPrefix = OutputPrefix + ".autosome15",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome16 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome16,
            AllCallShards = FilterMethylationShard.AllCallsAutosome16,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "16",
            OutputPrefix = OutputPrefix + ".autosome16",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome17 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome17,
            AllCallShards = FilterMethylationShard.AllCallsAutosome17,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "17",
            OutputPrefix = OutputPrefix + ".autosome17",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome18 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome18,
            AllCallShards = FilterMethylationShard.AllCallsAutosome18,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "18",
            OutputPrefix = OutputPrefix + ".autosome18",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome19 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome19,
            AllCallShards = FilterMethylationShard.AllCallsAutosome19,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "19",
            OutputPrefix = OutputPrefix + ".autosome19",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome20 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome20,
            AllCallShards = FilterMethylationShard.AllCallsAutosome20,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "20",
            OutputPrefix = OutputPrefix + ".autosome20",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome21 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome21,
            AllCallShards = FilterMethylationShard.AllCallsAutosome21,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "21",
            OutputPrefix = OutputPrefix + ".autosome21",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call MergeMethylationChromosome as MergeMethylationAutosome22 {
        input:
            FilteredCallShards = FilterMethylationShard.FilteredCallsAutosome22,
            AllCallShards = FilterMethylationShard.AllCallsAutosome22,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            Chromosome = AutosomePrefix + "22",
            OutputPrefix = OutputPrefix + ".autosome22",
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call AggregateMethylationChromosomes {
        input:
            FilteredCallsByChromosome = [MergeMethylationAutosome01.FilteredCalls, MergeMethylationAutosome02.FilteredCalls, MergeMethylationAutosome03.FilteredCalls, MergeMethylationAutosome04.FilteredCalls, MergeMethylationAutosome05.FilteredCalls, MergeMethylationAutosome06.FilteredCalls, MergeMethylationAutosome07.FilteredCalls, MergeMethylationAutosome08.FilteredCalls, MergeMethylationAutosome09.FilteredCalls, MergeMethylationAutosome10.FilteredCalls, MergeMethylationAutosome11.FilteredCalls, MergeMethylationAutosome12.FilteredCalls, MergeMethylationAutosome13.FilteredCalls, MergeMethylationAutosome14.FilteredCalls, MergeMethylationAutosome15.FilteredCalls, MergeMethylationAutosome16.FilteredCalls, MergeMethylationAutosome17.FilteredCalls, MergeMethylationAutosome18.FilteredCalls, MergeMethylationAutosome19.FilteredCalls, MergeMethylationAutosome20.FilteredCalls, MergeMethylationAutosome21.FilteredCalls, MergeMethylationAutosome22.FilteredCalls],
            SiteQCByChromosome = [MergeMethylationAutosome01.SiteQC, MergeMethylationAutosome02.SiteQC, MergeMethylationAutosome03.SiteQC, MergeMethylationAutosome04.SiteQC, MergeMethylationAutosome05.SiteQC, MergeMethylationAutosome06.SiteQC, MergeMethylationAutosome07.SiteQC, MergeMethylationAutosome08.SiteQC, MergeMethylationAutosome09.SiteQC, MergeMethylationAutosome10.SiteQC, MergeMethylationAutosome11.SiteQC, MergeMethylationAutosome12.SiteQC, MergeMethylationAutosome13.SiteQC, MergeMethylationAutosome14.SiteQC, MergeMethylationAutosome15.SiteQC, MergeMethylationAutosome16.SiteQC, MergeMethylationAutosome17.SiteQC, MergeMethylationAutosome18.SiteQC, MergeMethylationAutosome19.SiteQC, MergeMethylationAutosome20.SiteQC, MergeMethylationAutosome21.SiteQC, MergeMethylationAutosome22.SiteQC],
            SiteMetadataByChromosome = [MergeMethylationAutosome01.SiteMetadata, MergeMethylationAutosome02.SiteMetadata, MergeMethylationAutosome03.SiteMetadata, MergeMethylationAutosome04.SiteMetadata, MergeMethylationAutosome05.SiteMetadata, MergeMethylationAutosome06.SiteMetadata, MergeMethylationAutosome07.SiteMetadata, MergeMethylationAutosome08.SiteMetadata, MergeMethylationAutosome09.SiteMetadata, MergeMethylationAutosome10.SiteMetadata, MergeMethylationAutosome11.SiteMetadata, MergeMethylationAutosome12.SiteMetadata, MergeMethylationAutosome13.SiteMetadata, MergeMethylationAutosome14.SiteMetadata, MergeMethylationAutosome15.SiteMetadata, MergeMethylationAutosome16.SiteMetadata, MergeMethylationAutosome17.SiteMetadata, MergeMethylationAutosome18.SiteMetadata, MergeMethylationAutosome19.SiteMetadata, MergeMethylationAutosome20.SiteMetadata, MergeMethylationAutosome21.SiteMetadata, MergeMethylationAutosome22.SiteMetadata],
            RawMethylationBedByChromosome = [MergeMethylationAutosome01.RawMethylationBed, MergeMethylationAutosome02.RawMethylationBed, MergeMethylationAutosome03.RawMethylationBed, MergeMethylationAutosome04.RawMethylationBed, MergeMethylationAutosome05.RawMethylationBed, MergeMethylationAutosome06.RawMethylationBed, MergeMethylationAutosome07.RawMethylationBed, MergeMethylationAutosome08.RawMethylationBed, MergeMethylationAutosome09.RawMethylationBed, MergeMethylationAutosome10.RawMethylationBed, MergeMethylationAutosome11.RawMethylationBed, MergeMethylationAutosome12.RawMethylationBed, MergeMethylationAutosome13.RawMethylationBed, MergeMethylationAutosome14.RawMethylationBed, MergeMethylationAutosome15.RawMethylationBed, MergeMethylationAutosome16.RawMethylationBed, MergeMethylationAutosome17.RawMethylationBed, MergeMethylationAutosome18.RawMethylationBed, MergeMethylationAutosome19.RawMethylationBed, MergeMethylationAutosome20.RawMethylationBed, MergeMethylationAutosome21.RawMethylationBed, MergeMethylationAutosome22.RawMethylationBed],
            IntMethylationBedByChromosome = [MergeMethylationAutosome01.IntMethylationBed, MergeMethylationAutosome02.IntMethylationBed, MergeMethylationAutosome03.IntMethylationBed, MergeMethylationAutosome04.IntMethylationBed, MergeMethylationAutosome05.IntMethylationBed, MergeMethylationAutosome06.IntMethylationBed, MergeMethylationAutosome07.IntMethylationBed, MergeMethylationAutosome08.IntMethylationBed, MergeMethylationAutosome09.IntMethylationBed, MergeMethylationAutosome10.IntMethylationBed, MergeMethylationAutosome11.IntMethylationBed, MergeMethylationAutosome12.IntMethylationBed, MergeMethylationAutosome13.IntMethylationBed, MergeMethylationAutosome14.IntMethylationBed, MergeMethylationAutosome15.IntMethylationBed, MergeMethylationAutosome16.IntMethylationBed, MergeMethylationAutosome17.IntMethylationBed, MergeMethylationAutosome18.IntMethylationBed, MergeMethylationAutosome19.IntMethylationBed, MergeMethylationAutosome20.IntMethylationBed, MergeMethylationAutosome21.IntMethylationBed, MergeMethylationAutosome22.IntMethylationBed],
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            OutputPrefix = OutputPrefix,
            MemoryGB = MergeMemoryGB,
            DiskGB = MergeDiskGB,
            NumThreads = NumThreads
    }

    call AnnotateMethylationSites {
        input:
            SiteMetadata = AggregateMethylationChromosomes.SiteMetadata,
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
            BedFile = AggregateMethylationChromosomes.IntMethylationBed,
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
        File FilteredCalls = AggregateMethylationChromosomes.FilteredCalls
        File SiteQC = AggregateMethylationChromosomes.SiteQC
        File SiteMetadata = AggregateMethylationChromosomes.SiteMetadata
        File SampleQC = AggregateMethylationChromosomes.SampleQC
        File FilterSummary = AggregateMethylationChromosomes.FilterSummary
        File FilterCountsPlot = AggregateMethylationChromosomes.FilterCountsPlot
        File FilterUpsetPlot = AggregateMethylationChromosomes.FilterUpsetPlot
        File RawMethylationBed = AggregateMethylationChromosomes.RawMethylationBed
        File IntMethylationBed = AggregateMethylationChromosomes.IntMethylationBed
        File PassingSiteAnnotations = AnnotateMethylationSites.PassingSiteAnnotations
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        Array[File] ShardSampleQC = FilterMethylationShard.SampleQC
    }
}
