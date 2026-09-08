version 1.0
import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge

# Final phenotype PCs and optional TensorQTL covariate preparation.

workflow PrepareMethylationQtlCovariates {
    input {
        String calculate_phenotypepcs__computep_cs_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String mergecovariates__merge_covariatesr_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        File IntMethylationBed
        File? AdditionalCovariates
        String OutputPrefix
        Int PcMemoryGB
        Int PcDiskGB
        Int NumThreads

    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
      calculate_phenotypepcs__computep_cs_image = calculate_phenotypepcs__computep_cs_image,

            BedFile = IntMethylationBed,
            OutputPrefix = OutputPrefix + ".methylation",
            OutputSuffix = ".INT",
            memory = PcMemoryGB,
            disk_space = PcDiskGB,
            num_threads = NumThreads
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
    }

    output {
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
    }
}
