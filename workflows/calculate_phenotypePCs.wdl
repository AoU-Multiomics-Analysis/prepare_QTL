version 1.0

task ComputePCs{
    input {
        File BedFile
        String OutputPrefix
        String OutputSuffix = ""
        Int memory
        Int disk_space
        Int num_threads
    }
    command <<<

    Rscript /opt/prepare_qtl/scripts/common/calculate_PCs.R \
        --bed_file ~{BedFile} \
        --output_prefix ~{OutputPrefix} \
        --output_suffix ~{OutputSuffix}
    >>>
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File PhenotypePCsTSV="~{OutputPrefix}_phenotype_PCs~{OutputSuffix}.tsv"
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
    }
    call ComputePCs{
        input:
            BedFile = BedFile,
            OutputPrefix = OutputPrefix,
            OutputSuffix = OutputSuffix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    output {
        File OutPhenotypePCs= ComputePCs.PhenotypePCsTSV
    }

}
