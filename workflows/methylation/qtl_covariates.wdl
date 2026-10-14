version 1.0
import "../common/ExportRegenie.wdl" as Regenie
import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge

# Final phenotype PCs and optional TensorQTL covariate preparation.

workflow PrepareMethylationQtlCovariates {
    input {
        Int RegeniePhenotypesPerFile = 100
        String calculate_phenotypepcs__computep_cs_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        String mergecovariates__merge_covariatesr_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        File RawMethylationBed
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

    call Regenie.ExportRegenie as ExportIntRegenie {
        input:
            Bed = IntMethylationBed,
            Covariates = select_first([MergeIntAdditionalCovariates.QtlCovariates, IntPhenotypePCs.OutPhenotypePCs]),
            CovariateFormat = if defined(AdditionalCovariates) then "tensor" else "sample",
            OutputPrefix = OutputPrefix + ".methylation.INT",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = PcDiskGB
    }

    call Regenie.ExportRegenie as ExportRawRegenie {
        input:
            Bed = RawMethylationBed,
            Covariates = AdditionalCovariates,
            CovariateFormat = "sample",
            OutputPrefix = OutputPrefix + ".methylation.raw",
            PhenotypesPerFile = RegeniePhenotypesPerFile,
            DiskGB = PcDiskGB
    }

    output {
        Array[File] IntRegeniePhenotypes = ExportIntRegenie.Phenotypes
        File? IntRegenieCovariates = ExportIntRegenie.CovariatesOut
        File IntRegenieTraitMap = ExportIntRegenie.TraitMap
        File IntRegenieBatches = ExportIntRegenie.Batches
        File IntRegenieExcludedSamples = ExportIntRegenie.ExcludedSamples
        File IntRegenieLog = ExportIntRegenie.Log
        Array[File] RawRegeniePhenotypes = ExportRawRegenie.Phenotypes
        File? RawRegenieCovariates = ExportRawRegenie.CovariatesOut
        File RawRegenieTraitMap = ExportRawRegenie.TraitMap
        File RawRegenieBatches = ExportRawRegenie.Batches
        File RawRegenieExcludedSamples = ExportRawRegenie.ExcludedSamples
        File RawRegenieLog = ExportRawRegenie.Log

        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
    }
}
