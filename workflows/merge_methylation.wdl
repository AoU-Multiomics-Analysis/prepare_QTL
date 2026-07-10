version 1.0
import "calculate_phenotypePCs.wdl" as ComputePCs
import "MergeCovariates.wdl" as CovariateMerge

# The global merge is deliberately outside the scatter. A per-sample site
# filter would use a different denominator in every task and would therefore
# not implement MinSampleFraction across the cohort.

task FilterMethylationShard {
    input {
        String SampleId
        File MethylationBed
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

        printf 'sample_id\tfile_path\n%s\t%s\n' \
            "~{SampleId}" "~{MethylationBed}" > localized_manifest.tsv

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
        # TSV with sample_id and file_path columns. The workflow parses each
        # file_path into a typed File so Cromwell localizes gs:// objects.
        File SampleManifest
        String OutputPrefix
        File? AdditionalCovariates
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations

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

    Array[Map[String, String]] manifest_rows = read_objects(SampleManifest)

    scatter (sample_index in range(length(manifest_rows))) {
        Map[String, String] manifest_row = manifest_rows[sample_index]
        String sample_id = manifest_row["sample_id"]
        File methylation_bed = manifest_row["file_path"]
        String sample_output_prefix = "~{OutputPrefix}.sample.~{sample_index}"

        call FilterMethylationShard {
            input:
                SampleId = sample_id,
                MethylationBed = methylation_bed,
                OutputPrefix = sample_output_prefix,
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
            TotalSamples = length(manifest_rows),
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
