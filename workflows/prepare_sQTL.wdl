version 1.0
import  "calculate_phenotypePCs.wdl" as ComputePCs
import  "MergeCovariates.wdl" as CovariateMerge

task PrepareSpliceData {
    input {
        File SampleList
        File SpliceData
        String OutputPrefix

        Int memory
        Int disk_space
        Int num_threads
    }
    command {
        Rscript /tmp/PrepareSpliceData.R \
            --SpliceData ${SpliceData} \
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
        File IntSplicingBed="${OutputPrefix}.splicing.INT.bed.gz"
        File ScaledSplicingBed="${OutputPrefix}.splicing.scaled.bed.gz"
        File RawSplicingBed="${OutputPrefix}.splicing.raw.bed.gz"
        #File PhenotypeGroups = "${OutputPrefix}.phenotype_groups.tsv"
    }
 }

workflow sQTLPrepareData  {
    input {
        File SampleList
        File SpliceData
        String OutputPrefix
        File? AdditionalCovariates

        Int memory
        Int disk_space
        Int num_threads
    }
    call PrepareSpliceData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            SampleList = SampleList,
            SpliceData = SpliceData,
            OutputPrefix = OutputPrefix
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = PrepareSpliceData.IntSplicingBed,
            OutputPrefix = OutputPrefix + ".splicing",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = PrepareSpliceData.ScaledSplicingBed,
            OutputPrefix = OutputPrefix + ".splicing",
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
                OutputPrefix = OutputPrefix + ".splicing",
                OutputSuffix = ".INT"
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".splicing",
                OutputSuffix = ".scaled"
        }
    }
    output {
        File IntBedFile = PrepareSpliceData.IntSplicingBed
        File ScaledBedFile = PrepareSpliceData.ScaledSplicingBed
        File RawBedFile = PrepareSpliceData.RawSplicingBed
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File ScaledPhenotypePCsOut = ScaledPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        #File PhenotypeGroups = PrepareSpliceData.PhenotypeGroups
    }

}
