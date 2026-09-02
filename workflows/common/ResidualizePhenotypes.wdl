version 1.0

task ResidualizePhenotypes {
    input {
        File InputBed
        File? Covariates
        String OutputFileName

        Int memory
        Int disk_space
        Int num_threads
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2
    }

    command <<<
        set -euo pipefail
        stage="residualize_phenotypes"
        printf 'stage=%s start_time=%s dimensions=pending outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputFileName}"
        Rscript /tmp/ResidualizePhenotypes.R \
            --InputBed ~{InputBed} \
            ~{if defined(Covariates) then "--Covariates " + select_first([Covariates]) else ""} \
            --OutputFile ~{OutputFileName}
        printf 'stage=%s completion_time=%s dimensions=complete outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputFileName}"
    >>>

    runtime {
        docker: DockerImage
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
        preemptible: preemptible_attempts
        maxRetries: max_retries
    }

    output {
        File ResidualizedBed = "~{OutputFileName}"
    }
}
