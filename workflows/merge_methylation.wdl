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
        fwrite(manifest[, "sample_id", drop = FALSE], "cohort_samples.tsv", sep = "\t")
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
        File CohortSamples = "cohort_samples.tsv"
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

task MergeMethylationChromosome {
    input {
        Array[File] AllCallShards
        Array[File] SampleQCShards
        File CohortSamples
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
        printf '%s\n' ~{sep=' ' AllCallShards} > all_call_shards.list
        printf '%s\n' ~{sep=' ' SampleQCShards} > sample_qc_shards.list

        Rscript /tmp/MergeMethylationCohort.R \
            --AllCallList all_call_shards.list \
            --SampleQcList sample_qc_shards.list \
            --CohortSamples "~{CohortSamples}" \
            --TotalSamples ~{TotalSamples} \
            --Chromosome "~{Chromosome}" \
            --OutputPrefix "~{OutputPrefix}" \
            --MinSampleFraction ~{MinSampleFraction} \
            --MinSamples ~{MinSamples} \
            --MinMethylationMAD ~{MinMethylationMAD} \
            --ValueColumn "~{ValueColumn}" \
            --ValueMultiplier ~{ValueMultiplier} \
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

        concat_chromosome_tables() {
            local list_path="$1"
            local output_path="$2"
            local label="$3"
            local expected_header=""
            local input_path
            local current_header

            while IFS= read -r input_path; do
                [ -n "$input_path" ] || continue
                current_header=$(zgrep -m 1 '^' "$input_path")
                if [ -z "$expected_header" ]; then
                    expected_header="$current_header"
                elif [ "$current_header" != "$expected_header" ]; then
                    echo "${label} chromosome files do not have identical headers: ${input_path}" >&2
                    exit 1
                fi
            done < "$list_path"

            if [ -z "$expected_header" ]; then
                echo "${label} chromosome file list is empty" >&2
                exit 1
            fi

            {
                local first_file=1
                while IFS= read -r input_path; do
                    [ -n "$input_path" ] || continue
                    if [ "$first_file" -eq 1 ]; then
                        bgzip -c -d -@ ~{NumThreads} "$input_path"
                        first_file=0
                    else
                        bgzip -c -d -@ ~{NumThreads} "$input_path" | tail -n +2
                    fi
                done < "$list_path"
            } | bgzip -c -@ ~{NumThreads} > "$output_path"
        }

        # Inputs follow the autosome scatter order, and each chromosome task sorts
        # its own rows, so concatenation preserves genomic order without a second
        # whole-cohort sort.
        concat_chromosome_tables filtered_calls_by_chromosome.list \
            "~{OutputPrefix}.methylation.filtered.long.tsv.gz" "Filtered-call"
        concat_chromosome_tables site_qc_by_chromosome.list \
            "~{OutputPrefix}.methylation.site_qc.tsv.gz" "Site-QC"
        concat_chromosome_tables site_metadata_by_chromosome.list \
            "~{OutputPrefix}.methylation.site_metadata.tsv.gz" "Site-metadata"
        concat_chromosome_tables raw_beds_by_chromosome.list \
            "~{OutputPrefix}.methylation.raw.bed.gz" "Raw BED"
        concat_chromosome_tables int_beds_by_chromosome.list \
            "~{OutputPrefix}.methylation.INT.bed.gz" "INT BED"

        Rscript /tmp/AggregateMethylationChromosomes.R \
            --SiteMetadata "~{OutputPrefix}.methylation.site_metadata.tsv.gz" \
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

    Array[String] AutosomeNames = [
        AutosomePrefix + "1", AutosomePrefix + "2", AutosomePrefix + "3", AutosomePrefix + "4",
        AutosomePrefix + "5", AutosomePrefix + "6", AutosomePrefix + "7", AutosomePrefix + "8",
        AutosomePrefix + "9", AutosomePrefix + "10", AutosomePrefix + "11", AutosomePrefix + "12",
        AutosomePrefix + "13", AutosomePrefix + "14", AutosomePrefix + "15", AutosomePrefix + "16",
        AutosomePrefix + "17", AutosomePrefix + "18", AutosomePrefix + "19", AutosomePrefix + "20",
        AutosomePrefix + "21", AutosomePrefix + "22"
    ]
    Array[String] AutosomeOutputSuffixes = [
        "autosome01", "autosome02", "autosome03", "autosome04",
        "autosome05", "autosome06", "autosome07", "autosome08",
        "autosome09", "autosome10", "autosome11", "autosome12",
        "autosome13", "autosome14", "autosome15", "autosome16",
        "autosome17", "autosome18", "autosome19", "autosome20",
        "autosome21", "autosome22"
    ]
    Array[Array[File]] AllCallShardsByAutosome = transpose(FilterMethylationShard.AllCallsByAutosome)

    scatter (autosome_index in range(length(AutosomeNames))) {
        call MergeMethylationChromosome as MergeMethylationAutosome {
            input:
                AllCallShards = AllCallShardsByAutosome[autosome_index],
                SampleQCShards = FilterMethylationShard.SampleQC,
                CohortSamples = ShardMethylationManifest.CohortSamples,
                TotalSamples = ShardMethylationManifest.TotalSamples,
                Chromosome = AutosomeNames[autosome_index],
                OutputPrefix = OutputPrefix + "." + AutosomeOutputSuffixes[autosome_index],
                MinSampleFraction = MinSampleFraction,
                MinSamples = MinSamples,
                MinMethylationMAD = MinMethylationMAD,
                ValueColumn = ValueColumn,
                ValueMultiplier = ValueMultiplier,
                MemoryGB = MergeMemoryGB,
                DiskGB = MergeDiskGB,
                NumThreads = NumThreads
        }
    }

    call AggregateMethylationChromosomes {
        input:
            FilteredCallsByChromosome = MergeMethylationAutosome.FilteredCalls,
            SiteQCByChromosome = MergeMethylationAutosome.SiteQC,
            SiteMetadataByChromosome = MergeMethylationAutosome.SiteMetadata,
            RawMethylationBedByChromosome = MergeMethylationAutosome.RawMethylationBed,
            IntMethylationBedByChromosome = MergeMethylationAutosome.IntMethylationBed,
            SampleQcShards = FilterMethylationShard.SampleQC,
            TotalSamples = ShardMethylationManifest.TotalSamples,
            OutputPrefix = OutputPrefix,
            MemoryGB = AggregateMemoryGB,
            DiskGB = AggregateDiskGB,
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
