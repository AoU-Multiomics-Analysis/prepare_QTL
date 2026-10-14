version 1.0
import "../common/ExportRegenie.wdl" as Regenie

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
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
        docker: docker_image
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
        Int RegeniePhenotypesPerFile = 100
        String calculate_phenotypepcs__computep_cs_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String mergecovariates__merge_covariatesr_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String prepare_mqtl__prepare_methylation_data_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
            num_threads = num_threads,
            docker_image = prepare_mqtl__prepare_methylation_data_image
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,

            BedFile = PrepareMethylationData.IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
  }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,

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
      mergecovariates__merge_covariatesr_image = mergecovariates__merge_covariatesr_image,

                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT"
  }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
      mergecovariates__merge_covariatesr_image = mergecovariates__merge_covariatesr_image,

                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".scaled"
  }
    }

    call Regenie.ExportRegenie as ExportIntRegenie {
        input:
            Bed = PrepareMethylationData.IntMethylationBed,
            Covariates = select_first([MergeIntAdditionalCovariates.QtlCovariates, IntPhenotypePCs.OutPhenotypePCs]),
            CovariateFormat = if defined(AdditionalCovariates) then "tensor" else "sample",
            OutputPrefix = OutputPrefix + ".methylation.INT",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = disk_space
    }

    call Regenie.ExportRegenie as ExportScaledRegenie {
        input:
            Bed = PrepareMethylationData.ScaledMethylationBed,
            Covariates = select_first([MergeScaledAdditionalCovariates.QtlCovariates, ScaledPhenotypePCs.OutPhenotypePCs]),
            CovariateFormat = if defined(AdditionalCovariates) then "tensor" else "sample",
            OutputPrefix = OutputPrefix + ".methylation.scaled",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = disk_space
    }

    call Regenie.ExportRegenie as ExportRawRegenie {
        input:
            Bed = PrepareMethylationData.RawMethylationBed,
            Covariates = AdditionalCovariates,
            CovariateFormat = "sample",
            OutputPrefix = OutputPrefix + ".methylation.raw",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = disk_space
    }

    output {
        Array[File] IntRegeniePhenotypes = ExportIntRegenie.Phenotypes
        File? IntRegenieCovariates = ExportIntRegenie.CovariatesOut
        File IntRegenieTraitMap = ExportIntRegenie.TraitMap
        File IntRegenieBatches = ExportIntRegenie.Batches
        File IntRegenieExcludedSamples = ExportIntRegenie.ExcludedSamples
        File IntRegenieLog = ExportIntRegenie.Log
        Array[File] ScaledRegeniePhenotypes = ExportScaledRegenie.Phenotypes
        File? ScaledRegenieCovariates = ExportScaledRegenie.CovariatesOut
        File ScaledRegenieTraitMap = ExportScaledRegenie.TraitMap
        File ScaledRegenieBatches = ExportScaledRegenie.Batches
        File ScaledRegenieExcludedSamples = ExportScaledRegenie.ExcludedSamples
        File ScaledRegenieLog = ExportScaledRegenie.Log
        Array[File] RawRegeniePhenotypes = ExportRawRegenie.Phenotypes
        File? RawRegenieCovariates = ExportRawRegenie.CovariatesOut
        File RawRegenieTraitMap = ExportRawRegenie.TraitMap
        File RawRegenieBatches = ExportRawRegenie.Batches
        File RawRegenieExcludedSamples = ExportRawRegenie.ExcludedSamples
        File RawRegenieLog = ExportRawRegenie.Log

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
