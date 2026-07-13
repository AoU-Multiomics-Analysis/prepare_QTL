version 1.0 

workflow ConvertPlink {
    input {
        File vcf_file 
        Int new_id_max_allele_len
        String output_prefix
    }
    
    call plink2 {
        input:
            vcf_file = vcf_file,
            new_id_max_allele_len = new_id_max_allele_len,
            output_prefix = output_prefix
    }
    
    output  {
        File pgen = plink2.pgen 
        File pvar = plink2.pvar 
        File psam = plink2.psam
    }
}

task plink2 {
    input {
        File vcf_file
        String output_prefix
        Int new_id_max_allele_len
    }

    command <<<
        set -e

       # mkdir -p plink_output

        plink2 --vcf "~{vcf_file}" \
        --make-pgen \
        --out "~{output_prefix}" \
        --set-all-var-ids @:#_\$r_\$a \
        --new-id-max-allele-len "~{new_id_max_allele_len}" \
        --output-chr chrM \
        --chr 1-22
    >>>

    runtime {
        docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        memory: "96G"
        cpu: 4
        disks: "local-disk 400 SSD"
    }

    output {
        File pgen  = "~{output_prefix}.pgen"
        File pvar  = "~{output_prefix}.pvar"
        File psam  = "~{output_prefix}.psam"

    }
}



