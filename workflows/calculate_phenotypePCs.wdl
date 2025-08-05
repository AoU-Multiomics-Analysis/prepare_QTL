version 1.0 
workflow PhenotypePCs {
    input {
        File BedFile 
        String OutputPrefix
        Int memory 
        Int disk_space 
        Int num_threads
    }
    call ComputePCs{
        input:
            BedFile = BedFile,
            OutputPrefix = OutputPrefix,
            memory = memory, 
            disk_space = disk_space, 
            num_threads = num_threads
    } 

    output {
        File OutPhenotypePCs="${OutputPrefix}_phenotype_PCs.tsv"
    }

}

task ComputePCs{
    input {    
        File BedFile 
        String OutputPrefix 

        Int memory
        Int disk_space
        Int num_threads  
    }
    command <<<

    Rscript /tmp/compute_PCS.R \
        --expression_bed ${BedFile} \
        --output_prefix ${OutputPrefix}

    >>>
    runtime {
        docker: "AoU-Multiomics-Analysis/prepare_QTL:main"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File OutPhenotypePCs="${OutputPrefix}_phenotype_PCs.tsv"
    }
}

