version 1.0


workflow PrepareGenotypes{
    input {
        File InputVCF 
        File SampleList
        Int AlleleCountThreshold
        Int max_allele_len 
        String OutputPrefix

    }
    call SubsetVCF {
        input: 
            InputVCF = InputVCF,
            SampleList = SampleList, 
            OutputPrefix = OutputPrefix

    }
    
    call FilterVCF {
        input:
            VCF = SubsetVCF.OutSubsetVCF
            AlleleCountThreshold = AlleleCountThreshold,
            OutputPrefix = OutputPrefix
    }
    
    call plink2 {
        input:
            vcf_file = FilterVCF.FilteredVCF,
            output_prefix = OutputPrefix,
            new_id_max_allele_len = max_allele_len,
    }
    
    call ComputeGenotypePCs {
        vcf_file = ,
        output_prefix = OutputPrefix, 
        

    }
}


task SubsetVCF {
    input {
        File  InputVCF
        File SampleList
        String OutputPrefix
    }
    command <<< 
    bcftools view \
        -S ~{SampleList } \
        ~{InputVCF} \
        --Oz -output  ~{OutputPrefix}.vcf.gz 
    >>> 
    
    output {
        File OutSubsetVCF = "~{OutputPrefix}.vcf.bgz"
    }
    
    runtime {
        docker: ""
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }
}

task FilterVCF {
    input {
        File VCF
        Int AlleleCountThreshold
        String OutputPrefix
    }
    command <<<
    bcftools +fill-tags ~{VCF}  -- -t AC,AN,MAF,F_MISSING \
        | bcftools filter -i 'INFO/AC >= ~{AlleleCountThreshold}' \
        --Oz --output ~{OutputPrefix}.AC~{AlleleCountThreshold}.biallelic.vcf.gz 

    >>>
    output {
    File FilteredVCF = "~{OutputPrefix}.AC~{AlleleCountThreshold}.biallelic.vcf.gz"
    }

    runtime {
        docker: ""
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
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
        memory: "16G"
        cpu: 4
        disks: "local-disk 100 SSD"
    }

    output {
        Array[File] plink_outputs = glob("plink_output/*")
    }
}



task ComputeGenotypePCS {
    input {
        File vcf_file
        String output_prefix
        File genotype_rscript
    }

    command <<<
        set -e

        Rscript "~{genotype_rscript}" \
            --vcf_path "~{vcf_file}" \
            --prefix "~{output_prefix}"
        >>>
    
        runtime {
            docker: "quay.io/jonnguye/genotype_pcs:micromamba"
            memory: "8G"
            cpu: 2
        }
    
        output {
            File output_tsv = "~{output_prefix}_genetic_PCs.tsv"
        }
}


