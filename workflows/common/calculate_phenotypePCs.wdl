version 1.0

task ComputePCs{
    input {
        File BedFile
        String OutputPrefix
        String OutputSuffix = ""
        Int memory
        Int disk_space
        Int num_threads
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2
    }
    command <<<
        set -euo pipefail
        stage="compute_phenotype_pcs"
        printf 'stage=%s start_time=%s dimensions=pending outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.tsv,~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.all.tsv"
        Rscript /tmp/calculate_PCs.R \
            --bed_file ~{BedFile} \
            --output_prefix ~{OutputPrefix} \
            --output_suffix ~{OutputSuffix}
        printf 'stage=%s completion_time=%s dimensions=complete outputs=%s\n' \
            "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.tsv,~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.all.tsv"
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
        File PhenotypePCsTSV="~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.tsv"
        File PhenotypePCsAllTSV="~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.all.tsv"
    }
}



workflow PhenotypePCs {
    input {
        File BedFile
        String OutputPrefix
        String OutputSuffix = ""
        Int memory
        Int disk_space
        Int num_threads
        String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        Int preemptible_attempts = 2
        Int max_retries = 2
    }
    call ComputePCs{
        input:
            BedFile = BedFile,
            OutputPrefix = OutputPrefix,
            OutputSuffix = OutputSuffix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            DockerImage = DockerImage,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries
    }

    output {
        File OutPhenotypePCs= ComputePCs.PhenotypePCsTSV
        File OutPhenotypePCsAll = ComputePCs.PhenotypePCsAllTSV
    }

}
