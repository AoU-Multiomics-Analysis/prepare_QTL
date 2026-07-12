version 1.0
import "calculate_phenotypePCs.wdl" as ComputePCs
import "MergeCovariates.wdl" as CovariateMerge

# Cohort-level entry point for per-sample methylation outputs. In Terra, bind
# each Array[File] input to the corresponding output column over a sample-set
# (or other cohort) table entity.

task BuildMethylationCohortSamples {
    input {
        Array[File] SampleQCFiles
        String OutputPrefix
    }

    command <<<
        set -euo pipefail
        printf '%s\n' ~{sep=' ' SampleQCFiles} > sample_qc_files.list
        Rscript /tmp/BuildMethylationCohortSamples.R \
            --SampleQcList sample_qc_files.list \
            --OutputPrefix "~{OutputPrefix}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "4G"
        disks: "local-disk 20 HDD"
        cpu: 1
    }

    output {
        File CohortSamples = "~{OutputPrefix}.methylation.cohort_samples.tsv"
        File CohortSampleQC = "~{OutputPrefix}.methylation.cohort_sample_qc.tsv"
        Int TotalSamples = read_int("~{OutputPrefix}.methylation.total_samples.txt")
    }
}

task MergeMethylationChromosome {
    input {
        Array[File] AllCallShards
        File CohortSampleQC
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
        printf '%s\n' "~{CohortSampleQC}" > sample_qc_files.list

        Rscript /tmp/MergeMethylationCohort.R \
            --AllCallList all_call_shards.list \
            --SampleQcList sample_qc_files.list \
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
        File CohortSampleQC
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
        printf '%s\n' "~{CohortSampleQC}" > sample_qc_files.list

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
            --SampleQcList sample_qc_files.list \
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

workflow AggregateMethylationCohort {
    input {
        Array[File] AllCallsAutosome01
        Array[File] AllCallsAutosome02
        Array[File] AllCallsAutosome03
        Array[File] AllCallsAutosome04
        Array[File] AllCallsAutosome05
        Array[File] AllCallsAutosome06
        Array[File] AllCallsAutosome07
        Array[File] AllCallsAutosome08
        Array[File] AllCallsAutosome09
        Array[File] AllCallsAutosome10
        Array[File] AllCallsAutosome11
        Array[File] AllCallsAutosome12
        Array[File] AllCallsAutosome13
        Array[File] AllCallsAutosome14
        Array[File] AllCallsAutosome15
        Array[File] AllCallsAutosome16
        Array[File] AllCallsAutosome17
        Array[File] AllCallsAutosome18
        Array[File] AllCallsAutosome19
        Array[File] AllCallsAutosome20
        Array[File] AllCallsAutosome21
        Array[File] AllCallsAutosome22
        Array[File] SampleQCFiles
        String OutputPrefix
        File? AdditionalCovariates
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations

        Float MinSampleFraction = 0.95
        Int MinSamples = 0
        Float MinMethylationMAD = 0.003
        String AutosomePrefix = "chr"
        Int PromoterWindow = 2000
        String ValueColumn = "mod_score"
        Float ValueMultiplier = 0.01
        Int MergeMemoryGB = 128
        Int MergeDiskGB = 500
        Int AggregateMemoryGB = 64
        Int AggregateDiskGB = 1000
        Int AnnotationMemoryGB = 64
        Int AnnotationDiskGB = 100
        Int NumThreads = 1
    }

    call BuildMethylationCohortSamples {
        input:
            SampleQCFiles = SampleQCFiles,
            OutputPrefix = OutputPrefix
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
    Array[Array[File]] AllCallFilesByAutosome = [
        AllCallsAutosome01, AllCallsAutosome02, AllCallsAutosome03, AllCallsAutosome04,
        AllCallsAutosome05, AllCallsAutosome06, AllCallsAutosome07, AllCallsAutosome08,
        AllCallsAutosome09, AllCallsAutosome10, AllCallsAutosome11, AllCallsAutosome12,
        AllCallsAutosome13, AllCallsAutosome14, AllCallsAutosome15, AllCallsAutosome16,
        AllCallsAutosome17, AllCallsAutosome18, AllCallsAutosome19, AllCallsAutosome20,
        AllCallsAutosome21, AllCallsAutosome22
    ]

    scatter (autosome_index in range(length(AutosomeNames))) {
        call MergeMethylationChromosome as MergeMethylationAutosome {
            input:
                AllCallShards = AllCallFilesByAutosome[autosome_index],
                CohortSampleQC = BuildMethylationCohortSamples.CohortSampleQC,
                CohortSamples = BuildMethylationCohortSamples.CohortSamples,
                TotalSamples = BuildMethylationCohortSamples.TotalSamples,
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
            CohortSampleQC = BuildMethylationCohortSamples.CohortSampleQC,
            TotalSamples = BuildMethylationCohortSamples.TotalSamples,
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
        File CohortSamples = BuildMethylationCohortSamples.CohortSamples
        Int TotalSamples = BuildMethylationCohortSamples.TotalSamples
    }
}
