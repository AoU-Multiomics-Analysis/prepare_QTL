version 1.0
import  "calculate_phenotypePCs.wdl" as ComputePCs
import  "MergeCovariates.wdl" as CovariateMerge
import  "ResidualizePhenotypes.wdl" as Residualize

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
        Rscript /opt/prepare_qtl/scripts/splicing/PrepareSpliceData.R \
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

    if (ResidualizeNormalizedInputs) {
        call Residualize.ResidualizePhenotypes as ResidualizeIntPhenotypes {
            input:
                InputBed = PrepareSpliceData.IntSplicingBed,
                Covariates = MergeIntAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".splicing.INT.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = PrepareSpliceData.ScaledSplicingBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".splicing.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
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
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
        #File PhenotypeGroups = PrepareSpliceData.PhenotypeGroups
    }

}
