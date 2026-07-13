version 1.0
import "../common/calculate_phenotypePCs.wdl" as ComputePCs

# Phenotype-PC adjustment, local CpG correlation, and connectivity filtering.

task BuildMethylationCorrelationCovariates {
    input {
        File PhenotypePCs
        File? AdditionalCovariates
        String OutputPrefix
    }

    command <<<
        set -euo pipefail
        Rscript /tmp/BuildMethylationCorrelationCovariates.R \
            --PhenotypePCs "~{PhenotypePCs}" \
            ~{if defined(AdditionalCovariates) then "--AdditionalCovariates '" + select_first([AdditionalCovariates]) + "'" else ""} \
            --OutputFile "~{OutputPrefix}.methylation.correlation_covariates.tsv"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "16G"
        disks: "local-disk 50 HDD"
        cpu: 1
    }

    output {
        File CorrelationCovariates = "~{OutputPrefix}.methylation.correlation_covariates.tsv"
    }
}

task AnalyzeMethylationCpGCorrelation {
    input {
        File IntMethylationBed
        File Covariates
        String OutputPrefix
        Int WindowBP
        Float MinAbsCorrelation
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        Rscript /tmp/AnalyzeMethylationCpGCorrelation.R \
            --InputBed "~{IntMethylationBed}" \
            --Covariates "~{Covariates}" \
            --OutputPrefix "~{OutputPrefix}" \
            --WindowBP ~{WindowBP} \
            --MinAbsCorrelation ~{MinAbsCorrelation}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File CorrelationClusters = "~{OutputPrefix}.methylation.cpg_correlation_clusters.tsv.gz"
        File RepresentativeCpGs = "~{OutputPrefix}.methylation.cpg_correlation_representatives.tsv.gz"
        File CorrelationSummary = "~{OutputPrefix}.methylation.cpg_correlation_summary.tsv"
        File ClusterSizePlot = "~{OutputPrefix}.methylation.cpg_correlation_cluster_sizes.png"
        File MeanDistancePlot = "~{OutputPrefix}.methylation.cpg_correlation_mean_distances.png"
        File MaxDistancePlot = "~{OutputPrefix}.methylation.cpg_correlation_max_distances.png"
        File SpanPlot = "~{OutputPrefix}.methylation.cpg_correlation_span_vs_size.png"
    }
}


task FinalizeMethylationConnectivity {
    input {
        Array[File] RepresentativeCpGsByChromosome
        File PreConnectivityFilteredCalls
        File PreConnectivityRawMethylationBed
        File PreConnectivityIntMethylationBed
        File PreConnectivitySampleQC
        String OutputPrefix
        Float ConnectivityZThreshold
        Int MemoryGB
        Int DiskGB
    }

    command <<<
        set -euo pipefail
        printf '%s\n' "~{PreConnectivityIntMethylationBed}" > int_beds_by_chromosome.list
        printf '%s\n' ~{sep=' ' RepresentativeCpGsByChromosome} > representative_cpgs_by_chromosome.list

        Rscript /tmp/FinalizeMethylationConnectivity.R \
            --IntBedList int_beds_by_chromosome.list \
            --RepresentativeList representative_cpgs_by_chromosome.list \
            --FilteredCalls "~{PreConnectivityFilteredCalls}" \
            --RawBed "~{PreConnectivityRawMethylationBed}" \
            --IntBed "~{PreConnectivityIntMethylationBed}" \
            --SampleQC "~{PreConnectivitySampleQC}" \
            --OutputPrefix "~{OutputPrefix}" \
            --ConnectivityZThreshold ~{ConnectivityZThreshold}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File FilteredCalls = "~{OutputPrefix}.methylation.filtered.long.tsv.gz"
        File SampleQC = "~{OutputPrefix}.methylation.sample_qc.tsv"
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
        File ConnectivityOutliers = "~{OutputPrefix}.methylation.connectivity_outliers.tsv"
        File ConnectivitySummary = "~{OutputPrefix}.methylation.connectivity_summary.tsv"
        File ConnectivityRepresentativeCpGs = "~{OutputPrefix}.methylation.connectivity_representative_cpgs.tsv.gz"
    }
}


workflow RefineMethylationConnectivity {
    input {
        Array[File] IntMethylationBedsByChromosome
        Array[String] ChromosomeOutputSuffixes
        File PreConnectivityFilteredCalls
        File PreConnectivityRawMethylationBed
        File PreConnectivityIntMethylationBed
        File PreConnectivitySampleQC
        File? AdditionalCovariates
        String OutputPrefix
        Int PcMemoryGB = 128
        Int PcDiskGB = 500
        Int NumThreads = 1
        Int CorrelationWindowBP = 1000
        Float CorrelationMinAbsCorrelation = 0.95
        Int CorrelationMemoryGB = 64
        Int CorrelationDiskGB = 250
        Float ConnectivityZThreshold = -3.0
        Int ConnectivityMemoryGB = 64
        Int ConnectivityDiskGB = 1000
    }

    call ComputePCs.PhenotypePCs as PreliminaryIntPhenotypePCs {
        input:
            BedFile = PreConnectivityIntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation.pre_connectivity",
            OutputSuffix = ".INT",
            memory = PcMemoryGB,
            disk_space = PcDiskGB,
            num_threads = NumThreads
    }

    call BuildMethylationCorrelationCovariates {
        input:
            PhenotypePCs = PreliminaryIntPhenotypePCs.OutPhenotypePCs,
            AdditionalCovariates = AdditionalCovariates,
            OutputPrefix = OutputPrefix
    }

    scatter (chromosome_index in range(length(IntMethylationBedsByChromosome))) {
        call AnalyzeMethylationCpGCorrelation as AnalyzeMethylationChromosomeCorrelation {
            input:
                IntMethylationBed = IntMethylationBedsByChromosome[chromosome_index],
                Covariates = BuildMethylationCorrelationCovariates.CorrelationCovariates,
                OutputPrefix = OutputPrefix + "." + ChromosomeOutputSuffixes[chromosome_index],
                WindowBP = CorrelationWindowBP,
                MinAbsCorrelation = CorrelationMinAbsCorrelation,
                MemoryGB = CorrelationMemoryGB,
                DiskGB = CorrelationDiskGB
        }
    }

    call FinalizeMethylationConnectivity {
        input:
            RepresentativeCpGsByChromosome = AnalyzeMethylationChromosomeCorrelation.RepresentativeCpGs,
            PreConnectivityFilteredCalls = PreConnectivityFilteredCalls,
            PreConnectivityRawMethylationBed = PreConnectivityRawMethylationBed,
            PreConnectivityIntMethylationBed = PreConnectivityIntMethylationBed,
            PreConnectivitySampleQC = PreConnectivitySampleQC,
            OutputPrefix = OutputPrefix,
            ConnectivityZThreshold = ConnectivityZThreshold,
            MemoryGB = ConnectivityMemoryGB,
            DiskGB = ConnectivityDiskGB
    }

    output {
        File FilteredCalls = FinalizeMethylationConnectivity.FilteredCalls
        File SampleQC = FinalizeMethylationConnectivity.SampleQC
        File RawMethylationBed = FinalizeMethylationConnectivity.RawMethylationBed
        File IntMethylationBed = FinalizeMethylationConnectivity.IntMethylationBed
        File ConnectivityOutliers = FinalizeMethylationConnectivity.ConnectivityOutliers
        File ConnectivitySummary = FinalizeMethylationConnectivity.ConnectivitySummary
        File ConnectivityRepresentativeCpGs = FinalizeMethylationConnectivity.ConnectivityRepresentativeCpGs
        Array[File] CorrelationClustersByChromosome = AnalyzeMethylationChromosomeCorrelation.CorrelationClusters
        Array[File] CorrelationSummariesByChromosome = AnalyzeMethylationChromosomeCorrelation.CorrelationSummary
    }
}
