version 1.0
import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge

# Final phenotype PCs and optional TensorQTL covariate preparation.

workflow PrepareMethylationQtlCovariates {
    input {
        File IntMethylationBed
        File? AdditionalCovariates
        String OutputPrefix
        Int PcMemoryGB
        Int PcDiskGB
        Int NumThreads
        String methylation_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845"
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = PcMemoryGB,
            disk_space = PcDiskGB,
            num_threads = NumThreads,
            DockerImage = methylation_docker_image
    }

    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT",
                DockerImage = methylation_docker_image
        }
    }

    output {
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
    }
}
