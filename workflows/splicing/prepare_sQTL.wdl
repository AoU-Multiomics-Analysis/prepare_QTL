version 1.0
import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge
import "../common/ResidualizePhenotypes.wdl" as Residualize

task PrepareSpliceData {
    input {
        File SampleList
        File SpliceData
        String OutputPrefix

        Int memory
        Int disk_space
        Int num_threads
        String docker_image
    }
    command {
        Rscript /tmp/PrepareSpliceData.R \
            --SpliceData ${SpliceData} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}
        }

    runtime {
        docker: docker_image
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File IntSplicingBed="${OutputPrefix}.splicing.INT.bed.gz"
        File ScaledSplicingBed="${OutputPrefix}.splicing.scaled.bed.gz"
        File RawSplicingBed="${OutputPrefix}.splicing.raw.bed.gz"
        File IntConnectivityOutliers="${OutputPrefix}.splicing.INT.connectivity_outliers.tsv"
        File ScaledConnectivityOutliers="${OutputPrefix}.splicing.scaled.connectivity_outliers.tsv"
        #File PhenotypeGroups = "${OutputPrefix}.phenotype_groups.tsv"
    }
 }

workflow sQTLPrepareData  {
    input {
        File SampleList
        File SpliceData
        String OutputPrefix
        File? AdditionalCovariates
        Boolean ResidualizeNormalizedInputs = false

        Int memory
        Int disk_space
        Int num_threads
        String splicing_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845"
    }
    call PrepareSpliceData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            SampleList = SampleList,
            SpliceData = SpliceData,
            OutputPrefix = OutputPrefix,
            docker_image = splicing_docker_image
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = PrepareSpliceData.IntSplicingBed,
            OutputPrefix = OutputPrefix + ".splicing",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = splicing_docker_image
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = PrepareSpliceData.ScaledSplicingBed,
            OutputPrefix = OutputPrefix + ".splicing",
            OutputSuffix = ".scaled",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = splicing_docker_image
    }
    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".splicing",
                OutputSuffix = ".INT",
                DockerImage = splicing_docker_image
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".splicing",
                OutputSuffix = ".scaled",
                DockerImage = splicing_docker_image
        }
    }

    if (ResidualizeNormalizedInputs) {
        call Residualize.ResidualizePhenotypes as ResidualizeIntPhenotypes {
            input:
                InputBed = PrepareSpliceData.IntSplicingBed,
                Covariates = MergeIntAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".splicing.INT.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = splicing_docker_image
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = PrepareSpliceData.ScaledSplicingBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".splicing.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = splicing_docker_image
        }
    }

    output {
        File IntBedFile = PrepareSpliceData.IntSplicingBed
        File ScaledBedFile = PrepareSpliceData.ScaledSplicingBed
        File RawBedFile = PrepareSpliceData.RawSplicingBed
        File IntConnectivityOutliers = PrepareSpliceData.IntConnectivityOutliers
        File ScaledConnectivityOutliers = PrepareSpliceData.ScaledConnectivityOutliers
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File ScaledPhenotypePCsOut = ScaledPhenotypePCs.OutPhenotypePCs
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File ScaledPhenotypePCsAllOut = ScaledPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
        #File PhenotypeGroups = PrepareSpliceData.PhenotypeGroups
    }

}
