version 1.0

import  "calculate_phenotypePCs.wdl" as ComputePCs
import  "MergeCovariates.wdl" as CovariateMerge
import  "ResidualizePhenotypes.wdl" as Residualize



task eqtl_prepare_expression {
    input {
        File CountGCT
        File AnnotationGTF
        File SampleList
        String OutputPrefix


        Int memory
        Int disk_space
        Int num_threads

        }
    command {
        Rscript /opt/prepare_qtl/scripts/expression/PrepareExpression.R \
            --CountGCT ${CountGCT} \
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
        File IntExpressionBed="${OutputPrefix}.expression.INT.bed.gz"
        File ScaledExpressionBed="${OutputPrefix}.expression.scaled.bed.gz"
        File RawExpressionBed="${OutputPrefix}.expression.raw.bed.gz"
        File IntConnectivityOutliers="${OutputPrefix}.expression.INT.connectivity_outliers.tsv"
        File ScaledConnectivityOutliers="${OutputPrefix}.expression.scaled.connectivity_outliers.tsv"
    }
}

workflow eQTLPrepareData {
    input {
        String OutputPrefix
        File CountGCT
        File AnnotationGTF
        File SampleList
        File? AdditionalCovariates
        Boolean ResidualizeNormalizedInputs = false

        Int memory
        Int disk_space
        Int num_threads

            }
    call eqtl_prepare_expression {
        input:
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            CountGCT  = CountGCT,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = eqtl_prepare_expression.IntExpressionBed,
            OutputPrefix = OutputPrefix + ".expression",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = eqtl_prepare_expression.ScaledExpressionBed,
            OutputPrefix = OutputPrefix + ".expression",
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
                OutputPrefix = OutputPrefix + ".expression",
                OutputSuffix = ".INT"
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".expression",
                OutputSuffix = ".scaled"
        }
    }

    if (ResidualizeNormalizedInputs) {
        call Residualize.ResidualizePhenotypes as ResidualizeIntPhenotypes {
            input:
                InputBed = eqtl_prepare_expression.IntExpressionBed,
                Covariates = MergeIntAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".expression.INT.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = eqtl_prepare_expression.ScaledExpressionBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".expression.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
        }
    }

    output {
        File IntBedFile = eqtl_prepare_expression.IntExpressionBed
        File ScaledBedFile = eqtl_prepare_expression.ScaledExpressionBed
        File RawBedFile = eqtl_prepare_expression.RawExpressionBed
        File IntConnectivityOutliers = eqtl_prepare_expression.IntConnectivityOutliers
        File ScaledConnectivityOutliers = eqtl_prepare_expression.ScaledConnectivityOutliers
        File IntPhenotypePCsOut = IntPhenotypePCs.OutPhenotypePCs
        File ScaledPhenotypePCsOut = ScaledPhenotypePCs.OutPhenotypePCs
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
    }
}
