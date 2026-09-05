version 1.0
import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge
import "../common/ResidualizePhenotypes.wdl" as Residualize




task PrepareProteomicData {
    input {
        File AnnotationGTF
        File SampleList
        File ProteomicData
        String OutputPrefix

        Int memory
        Int disk_space
        Int num_threads
        String docker_image
    }
    command {
        Rscript /tmp/PrepareProteomics.R \
            --ProteomicData ${ProteomicData} \
            --AnnotationGTF ${AnnotationGTF} \
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
        File IntProteomicBed="${OutputPrefix}.protein.INT.bed.gz"
        File ScaledProteomicBed="${OutputPrefix}.protein.scaled.bed.gz"
        File RawProteomicBed="${OutputPrefix}.protein.raw.bed.gz"
        File IntConnectivityOutliers="${OutputPrefix}.protein.INT.connectivity_outliers.tsv"
        File ScaledConnectivityOutliers="${OutputPrefix}.protein.scaled.connectivity_outliers.tsv"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow pQTLPrepareData {
    input {
        Int memory
        Int disk_space
        Int num_threads
        File AnnotationGTF
        File SampleList
        File ProteomicData
        String OutputPrefix
        File? AdditionalCovariates
        Boolean ResidualizeNormalizedInputs = false
        String proteomics_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845"
    }
    call PrepareProteomicData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix,
            ProteomicData = ProteomicData,
            docker_image = proteomics_docker_image

    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = PrepareProteomicData.IntProteomicBed,
            OutputPrefix = OutputPrefix + ".protein",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = proteomics_docker_image
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = PrepareProteomicData.ScaledProteomicBed,
            OutputPrefix = OutputPrefix + ".protein",
            OutputSuffix = ".scaled",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = proteomics_docker_image
    }
    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".protein",
                OutputSuffix = ".INT",
                DockerImage = proteomics_docker_image
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".protein",
                OutputSuffix = ".scaled",
                DockerImage = proteomics_docker_image
        }
    }

    if (ResidualizeNormalizedInputs) {
        call Residualize.ResidualizePhenotypes as ResidualizeIntPhenotypes {
            input:
                InputBed = PrepareProteomicData.IntProteomicBed,
                Covariates = MergeIntAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".protein.INT.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = proteomics_docker_image
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = PrepareProteomicData.ScaledProteomicBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".protein.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = proteomics_docker_image
        }
    }

    output {
        File IntBedFile = PrepareProteomicData.IntProteomicBed
        File ScaledBedFile = PrepareProteomicData.ScaledProteomicBed
        File RawBedFile = PrepareProteomicData.RawProteomicBed
        File IntConnectivityOutliers = PrepareProteomicData.IntConnectivityOutliers
        File ScaledConnectivityOutliers = PrepareProteomicData.ScaledConnectivityOutliers
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File ScaledPhenotypePCsOut = ScaledPhenotypePCs.OutPhenotypePCs
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File ScaledPhenotypePCsAllOut = ScaledPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
    }
}
