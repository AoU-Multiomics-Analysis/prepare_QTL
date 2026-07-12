version 1.0
import  "calculate_phenotypePCs.wdl" as ComputePCs
import  "MergeCovariates.wdl" as CovariateMerge
import  "ResidualizePhenotypes.wdl" as Residualize




task PrepareProteomicData {
    input {
        File AnnotationGTF
        File SampleList
        File ProteomicData
        String OutputPrefix

        Int memory
        Int disk_space
        Int num_threads
    }
    command {
        Rscript /opt/prepare_qtl/scripts/proteomics/PrepareProteomics.R \
            --ProteomicData ${ProteomicData} \
            --AnnotationGTF ${AnnotationGTF} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}
        }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
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
    }
    call PrepareProteomicData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix,
            ProteomicData = ProteomicData

    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = PrepareProteomicData.IntProteomicBed,
            OutputPrefix = OutputPrefix + ".protein",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = PrepareProteomicData.ScaledProteomicBed,
            OutputPrefix = OutputPrefix + ".protein",
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
                OutputPrefix = OutputPrefix + ".protein",
                OutputSuffix = ".INT"
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".protein",
                OutputSuffix = ".scaled"
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
                num_threads = num_threads
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = PrepareProteomicData.ScaledProteomicBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".protein.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
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
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
    }
}
