version 1.0
import "../common/ExportRegenie.wdl" as Regenie
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
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
        Int RegeniePhenotypesPerFile = 100
        String calculate_phenotypepcs__computep_cs_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String mergecovariates__merge_covariatesr_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String prepare_pqtl__prepare_proteomic_data_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String residualizephenotypes__residualize_phenotypes_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
            ProteomicData = ProteomicData,
            docker_image = prepare_pqtl__prepare_proteomic_data_image

    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,

            BedFile = PrepareProteomicData.IntProteomicBed,
            OutputPrefix = OutputPrefix + ".protein",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
  }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,

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
      mergecovariates__merge_covariatesr_image = mergecovariates__merge_covariatesr_image,

                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".protein",
                OutputSuffix = ".INT"
  }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
      mergecovariates__merge_covariatesr_image = mergecovariates__merge_covariatesr_image,

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
                num_threads = num_threads,
                DockerImage = residualizephenotypes__residualize_phenotypes_image
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = PrepareProteomicData.ScaledProteomicBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".protein.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = residualizephenotypes__residualize_phenotypes_image
        }
    }

    call Regenie.ExportRegenie as ExportIntRegenie {
        input:
            Bed = PrepareProteomicData.IntProteomicBed,
            Covariates = select_first([MergeIntAdditionalCovariates.QtlCovariates, IntPhenotypePCs.OutPhenotypePCs]),
            CovariateFormat = if defined(AdditionalCovariates) then "tensor" else "sample",
            OutputPrefix = OutputPrefix + ".protein.INT",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = disk_space
    }

    call Regenie.ExportRegenie as ExportScaledRegenie {
        input:
            Bed = PrepareProteomicData.ScaledProteomicBed,
            Covariates = select_first([MergeScaledAdditionalCovariates.QtlCovariates, ScaledPhenotypePCs.OutPhenotypePCs]),
            CovariateFormat = if defined(AdditionalCovariates) then "tensor" else "sample",
            OutputPrefix = OutputPrefix + ".protein.scaled",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = disk_space
    }

    call Regenie.ExportRegenie as ExportRawRegenie {
        input:
            Bed = PrepareProteomicData.RawProteomicBed,
            Covariates = AdditionalCovariates,
            CovariateFormat = "sample",
            OutputPrefix = OutputPrefix + ".protein.raw",
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
