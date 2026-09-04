version 1.0

import "../common/calculate_phenotypePCs.wdl" as ComputePCs
import "../common/MergeCovariates.wdl" as CovariateMerge
import "../common/ResidualizePhenotypes.wdl" as Residualize



task eqtl_prepare_expression {
    input {
        File? CountGCT
        File? AnnotationGTF
        File? CpmBed
        File? Log2CpmBed
        File SampleList
        String OutputPrefix


        Int memory
        Int disk_space
        Int num_threads
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2

        }
    command {
        set -euo pipefail
        stage="eqtl_prepare_expression"
        printf 'stage=%s start_time=%s dimensions=pending outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}.expression.INT.bed.gz,~{OutputPrefix}.expression.scaled.bed.gz,~{OutputPrefix}.expression.raw.bed.gz,~{OutputPrefix}.expression.INT.connectivity_outliers.tsv,~{OutputPrefix}.expression.scaled.connectivity_outliers.tsv"
        Rscript /tmp/PrepareExpression.R \
            ~{if defined(CountGCT) then "--CountGCT \"" + select_first([CountGCT]) + "\"" else ""} \
            ~{if defined(AnnotationGTF) then "--AnnotationGTF \"" + select_first([AnnotationGTF]) + "\"" else ""} \
            ~{if defined(CpmBed) then "--CpmBed '" + sub(select_first([CpmBed]), "'", "'\"'\"'") + "'" else ""} \
            ~{if defined(Log2CpmBed) then "--Log2CpmBed \"" + select_first([Log2CpmBed]) + "\"" else ""} \
            --SampleList "~{SampleList}" \
            --OutputPrefix "~{OutputPrefix}"
        printf 'stage=%s completion_time=%s dimensions=complete outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}.expression.INT.bed.gz,~{OutputPrefix}.expression.scaled.bed.gz,~{OutputPrefix}.expression.raw.bed.gz,~{OutputPrefix}.expression.INT.connectivity_outliers.tsv,~{OutputPrefix}.expression.scaled.connectivity_outliers.tsv"

        }

    runtime {
        docker: DockerImage
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: preemptible_attempts
        maxRetries: max_retries
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
        File? CountGCT
        File? AnnotationGTF
        File? CpmBed
        File? Log2CpmBed
        File SampleList
        File? AdditionalCovariates
        Boolean ResidualizeNormalizedInputs = false

        Int memory
        Int disk_space
        Int num_threads
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2

            }
    call eqtl_prepare_expression {
        input:
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = DockerImage,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries,
            CountGCT  = CountGCT,
            AnnotationGTF = AnnotationGTF,
            CpmBed = CpmBed,
            Log2CpmBed = Log2CpmBed,
            SampleList = SampleList
    }

    call ComputePCs.PhenotypePCs as IntPhenotypePCs {
        input:
            BedFile = eqtl_prepare_expression.IntExpressionBed,
            OutputPrefix = OutputPrefix + ".expression",
            OutputSuffix = ".INT",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = DockerImage,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries
    }

    call ComputePCs.PhenotypePCs as ScaledPhenotypePCs {
        input:
            BedFile = eqtl_prepare_expression.ScaledExpressionBed,
            OutputPrefix = OutputPrefix + ".expression",
            OutputSuffix = ".scaled",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = DockerImage,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries
    }

    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeIntAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = IntPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".expression",
                OutputSuffix = ".INT",
                DockerImage = DockerImage,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries
        }

        call CovariateMerge.MergeCovariates as MergeScaledAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = ScaledPhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix + ".expression",
                OutputSuffix = ".scaled",
                DockerImage = DockerImage,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries
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
                num_threads = num_threads,
                DockerImage = DockerImage,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries
        }

        call Residualize.ResidualizePhenotypes as ResidualizeScaledPhenotypes {
            input:
                InputBed = eqtl_prepare_expression.ScaledExpressionBed,
                Covariates = MergeScaledAdditionalCovariates.QtlCovariates,
                OutputFileName = OutputPrefix + ".expression.scaled.residualized.bed.gz",
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads,
                DockerImage = DockerImage,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries
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
        File IntPhenotypePCsAllOut = IntPhenotypePCs.OutPhenotypePCsAll
        File ScaledPhenotypePCsAllOut = ScaledPhenotypePCs.OutPhenotypePCsAll
        File? IntQtlCovariates = MergeIntAdditionalCovariates.QtlCovariates
        File? ScaledQtlCovariates = MergeScaledAdditionalCovariates.QtlCovariates
        File? IntResidualizedBedFile = ResidualizeIntPhenotypes.ResidualizedBed
        File? ScaledResidualizedBedFile = ResidualizeScaledPhenotypes.ResidualizedBed
    }
}
