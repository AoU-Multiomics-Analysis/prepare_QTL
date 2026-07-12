version 1.0
import "AggregateMethylationCohort.wdl" as CohortAggregation

# Manifest/shard entry point. The global cohort reduction is delegated to
# AggregateMethylationCohort.wdl so it is shared with Terra-table processing.

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

        if [ ~{NumThreads} -lt 1 ]; then
            echo "NumThreads must be at least 1" >&2
            exit 1
        fi

        mkdir -p input_beds
        printf 'sample_id\tfile_path\n' > localized_manifest.tsv
        : > transfer_args.bin

        tail -n +2 "~{ShardManifest}" | while IFS=$'\t' read -r sample_id source_path; do
            [ -n "$sample_id" ] || continue
            local_path="input_beds/${sample_id}.combined.bed.gz"
            printf '%s\t%s\n' "$sample_id" "$local_path" >> localized_manifest.tsv
            printf '%s\0%s\0%s\0' "$sample_id" "$source_path" "$local_path" >> transfer_args.bin
        done

        # shellcheck disable=SC2016
        xargs -0 -n 3 -P ~{NumThreads} bash -c '
            sample_id="$1"
            source_path="$2"
            local_path="$3"
            if [[ "$source_path" == gs://* ]]; then
                gsutil -q cp "$source_path" "$local_path"
            else
                if [ ! -f "$source_path" ]; then
                    echo "Input BED file for ${sample_id} is not accessible inside the task: ${source_path}" >&2
                    exit 1
                fi
                cp "$source_path" "$local_path"
            fi
        ' _ < transfer_args.bin

        Rscript /tmp/FilterMethylationShard.R \
            --InputManifest localized_manifest.tsv \
            --OutputPrefix "~{OutputPrefix}" \
            --MinCoverage ~{MinCoverage} \
            --FilterChroms "~{FilterChroms}" \
            --FenceK ~{FenceK} \
            --AutosomePrefix "~{AutosomePrefix}" \
            --NumThreads ~{NumThreads}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: "~{NumThreads}"
    }

    output {
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
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
        Array[File] AllCallsByAutosome = [
            AllCallsAutosome01, AllCallsAutosome02, AllCallsAutosome03, AllCallsAutosome04,
            AllCallsAutosome05, AllCallsAutosome06, AllCallsAutosome07, AllCallsAutosome08,
            AllCallsAutosome09, AllCallsAutosome10, AllCallsAutosome11, AllCallsAutosome12,
            AllCallsAutosome13, AllCallsAutosome14, AllCallsAutosome15, AllCallsAutosome16,
            AllCallsAutosome17, AllCallsAutosome18, AllCallsAutosome19, AllCallsAutosome20,
            AllCallsAutosome21, AllCallsAutosome22
        ]
    }
}

workflow MergeMethylation {
    input {
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
        Int ShardMemoryGB = 64
        Int ShardDiskGB = 250
        Int ShardNumThreads = 4
        Int MergeMemoryGB = 128
        Int MergeDiskGB = 500
        Int AggregateMemoryGB = 64
        Int AggregateDiskGB = 1000
        Int AnnotationMemoryGB = 64
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
                NumThreads = ShardNumThreads
        }
    }

    Array[Array[File]] AllCallShardsByAutosome = transpose(FilterMethylationShard.AllCallsByAutosome)

    call CohortAggregation.AggregateMethylationCohort as CohortMerge {
        input:
            AllCallsAutosome01 = AllCallShardsByAutosome[0],
            AllCallsAutosome02 = AllCallShardsByAutosome[1],
            AllCallsAutosome03 = AllCallShardsByAutosome[2],
            AllCallsAutosome04 = AllCallShardsByAutosome[3],
            AllCallsAutosome05 = AllCallShardsByAutosome[4],
            AllCallsAutosome06 = AllCallShardsByAutosome[5],
            AllCallsAutosome07 = AllCallShardsByAutosome[6],
            AllCallsAutosome08 = AllCallShardsByAutosome[7],
            AllCallsAutosome09 = AllCallShardsByAutosome[8],
            AllCallsAutosome10 = AllCallShardsByAutosome[9],
            AllCallsAutosome11 = AllCallShardsByAutosome[10],
            AllCallsAutosome12 = AllCallShardsByAutosome[11],
            AllCallsAutosome13 = AllCallShardsByAutosome[12],
            AllCallsAutosome14 = AllCallShardsByAutosome[13],
            AllCallsAutosome15 = AllCallShardsByAutosome[14],
            AllCallsAutosome16 = AllCallShardsByAutosome[15],
            AllCallsAutosome17 = AllCallShardsByAutosome[16],
            AllCallsAutosome18 = AllCallShardsByAutosome[17],
            AllCallsAutosome19 = AllCallShardsByAutosome[18],
            AllCallsAutosome20 = AllCallShardsByAutosome[19],
            AllCallsAutosome21 = AllCallShardsByAutosome[20],
            AllCallsAutosome22 = AllCallShardsByAutosome[21],
            SampleQCFiles = FilterMethylationShard.SampleQC,
            OutputPrefix = OutputPrefix,
            AdditionalCovariates = AdditionalCovariates,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            MinSampleFraction = MinSampleFraction,
            MinSamples = MinSamples,
            MinMethylationMAD = MinMethylationMAD,
            AutosomePrefix = AutosomePrefix,
            PromoterWindow = PromoterWindow,
            ValueColumn = ValueColumn,
            ValueMultiplier = ValueMultiplier,
            MergeMemoryGB = MergeMemoryGB,
            MergeDiskGB = MergeDiskGB,
            AggregateMemoryGB = AggregateMemoryGB,
            AggregateDiskGB = AggregateDiskGB,
            AnnotationMemoryGB = AnnotationMemoryGB,
            AnnotationDiskGB = AnnotationDiskGB,
            NumThreads = NumThreads
    }

    output {
        File FilteredCalls = CohortMerge.FilteredCalls
        File SiteQC = CohortMerge.SiteQC
        File SiteMetadata = CohortMerge.SiteMetadata
        File SampleQC = CohortMerge.SampleQC
        File FilterSummary = CohortMerge.FilterSummary
        File FilterCountsPlot = CohortMerge.FilterCountsPlot
        File FilterUpsetPlot = CohortMerge.FilterUpsetPlot
        File RawMethylationBed = CohortMerge.RawMethylationBed
        File IntMethylationBed = CohortMerge.IntMethylationBed
        File PassingSiteAnnotations = CohortMerge.PassingSiteAnnotations
        File IntPhenotypePCsOut = CohortMerge.IntPhenotypePCsOut
        File? IntQtlCovariates = CohortMerge.IntQtlCovariates
        Array[File] ShardSampleQC = FilterMethylationShard.SampleQC
    }
}
