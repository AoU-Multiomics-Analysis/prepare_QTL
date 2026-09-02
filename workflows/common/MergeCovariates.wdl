version 1.0


workflow MergeCovariates {
    input {
        File GenotypePCs
        String OutputPrefix
        String OutputSuffix = ""
        File MolecularPCs
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2
    }

    call MergeCovariatesR {
        input:
            GenotypePCs = GenotypePCs,
            OutputPrefix = OutputPrefix,
            OutputSuffix = OutputSuffix,
            MolecularPCs = MolecularPCs,
            DockerImage = DockerImage,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries
    }

    output {
        File QtlCovariates = MergeCovariatesR.QtlCovariates

    }
}


task MergeCovariatesR {
    input {
        File GenotypePCs
        String OutputPrefix
        String OutputSuffix = ""
        File MolecularPCs
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2
    }

    command <<<
        set -euo pipefail
        stage="merge_covariates"
        printf 'stage=%s start_time=%s dimensions=pending outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}_QTL_covariates~{OutputSuffix}.tsv"
        Rscript /tmp/MergeCovariates.R \
            --GenotypePCs ~{GenotypePCs} \
            --MolecularPCs ~{MolecularPCs} \
            --OutputPrefix ~{OutputPrefix} \
            --OutputSuffix ~{OutputSuffix}
        printf 'stage=%s completion_time=%s dimensions=complete outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}_QTL_covariates~{OutputSuffix}.tsv"
        >>>

        runtime {
            docker: DockerImage
            memory: "96G"
            cpu: 1
            disks: "local-disk 100 SSD"
            preemptible: preemptible_attempts
            maxRetries: max_retries
        }

        output {
            File QtlCovariates = "~{OutputPrefix}_QTL_covariates~{OutputSuffix}.tsv"
        }



}
