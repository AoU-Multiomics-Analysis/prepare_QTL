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
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
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
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".methylation",
                OutputSuffix = ".INT"
        }
    }

    output {
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
    }
}
