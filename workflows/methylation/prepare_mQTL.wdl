version 1.0

import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge

task PrepareMethylationData {
    input {
        File MethylationBed
        File SampleList
        String OutputPrefix
        Float MissingnessThreshold = 5.0

        Int memory
        Int disk_space
        Int num_threads
    }

    command <<<
        set -euo pipefail
        echo "Starting methylation BED preparation for ~{MethylationBed}"
        echo "Using a strict feature missingness threshold of < ~{MissingnessThreshold}%"
        Rscript /tmp/PrepareMethylation.R \
            --MethylationBed "~{MethylationBed}" \
            --SampleList "~{SampleList}" \
            --OutputPrefix "~{OutputPrefix}" \
            --MissingnessThreshold ~{MissingnessThreshold}
        echo "Completed methylation BED preparation for ~{OutputPrefix}"
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File IntMethylationBed = "~{OutputPrefix}.methylation.INT.bed.gz"
        File ScaledMethylationBed = "~{OutputPrefix}.methylation.scaled.bed.gz"
        File RawMethylationBed = "~{OutputPrefix}.methylation.raw.bed.gz"
        File IntConnectivityOutliers = "~{OutputPrefix}.methylation.INT.connectivity_outliers.tsv"
        File ScaledConnectivityOutliers = "~{OutputPrefix}.methylation.scaled.connectivity_outliers.tsv"
    }
}

workflow mQTLPrepareData {
    input {
        File MethylationBed
        File SampleList
        String OutputPrefix
        File? AdditionalCovariates
        Float MissingnessThreshold = 5.0

        Int memory
        Int disk_space
        Int num_threads
    }

    call PrepareMethylationData {
        input:
            MethylationBed = MethylationBed,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix,
            MissingnessThreshold = MissingnessThreshold,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = PrepareMethylationData.IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = PrepareMethylationData.ScaledMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".scaled",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT"
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".scaled"
        }
    }

    output {
        File IntBedFile = PrepareMethylationData.IntMethylationBed
        File ScaledBedFile = PrepareMethylationData.ScaledMethylationBed
        File RawBedFile = PrepareMethylationData.RawMethylationBed
        File IntConnectivityOutliers = PrepareMethylationData.IntConnectivityOutliers
        File ScaledConnectivityOutliers = PrepareMethylationData.ScaledConnectivityOutliers
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File ScaledPhenotypePCsOut = ScaledPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
    }
}
