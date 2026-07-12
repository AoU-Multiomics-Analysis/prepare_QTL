version 1.0

task ResidualizePhenotypes {
    input {
        File InputBed
        File? Covariates
        String OutputFileName

        Int memory
        Int disk_space
        Int num_threads
    }

    command <<<
        Rscript /opt/prepare_qtl/scripts/common/ResidualizePhenotypes.R \
            --InputBed ~{InputBed} \
            ~{if defined(Covariates) then "--Covariates " + select_first([Covariates]) else ""} \
            --OutputFile ~{OutputFileName}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File ResidualizedBed = "~{OutputFileName}"
    }
}
