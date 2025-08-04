task PrepareProteomicData {
        File AnnotationGTF
        File SampleList 
        File ProteomicData 
        String OutputPrefix 
        
        Int memory 
        Int disk_space 
        Int num_threads 
    command {
        Rscript /tmp/PrepareProteomics.R \
            --ProteomicData ${ProteomicData} \
            --AnnotationGTF ${AnnotationGTF} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}
        }

    runtime {
        docker: "quay.io/jonnguye/modified_gtex_eqtl:1.1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File expression_bed="${prefix}.expression.bed.gz"
        File expression_bed_index="${prefix}.expression.bed.gz.tbi"
    }

    meta {
        author: "Francois Aguet"
    }
}

task compute_PCs{
    
        File expression_bed 
        File genotype_covariates 
        String prefix

        Int memory
        Int disk_space
        Int num_threads   
    command {

    Rscript /tmp/compute_PCS.R \
        --expression_bed ${expression_bed} \
        --genotype_covariates ${genotype_covariates} \
        --output_prefix ${prefix}

    }
    runtime {
        docker: "evinpadhi/prepare_qtl:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File phenotype_pcs="${prefix}_phenotype_PCs.tsv"
        File QTL_covariates="${prefix}_QTL_covariates.tsv"
    }

    meta {
        author: "Evin Padhi"
    }
} 

workflow prepare_QTL_data {
        String prefix
        Int memory 
        Int disk_space 
        Int num_threads 
        File genotype_covariates
        

        File tpm_gct
        File counts_gct
        File annotation_gtf
        File sample_participant_ids
        File vcf_chr_list
        File sample_list
        Float? tpm_threshold
        Int? count_threshold
        Float? sample_frac_threshold
        String? normalization_method
        String? flags
    call eqtl_prepare_expression {
        input:
            prefix = prefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            tpm_gct = tpm_gct,
            counts_gct = counts_gct,
            annotation_gtf = annotation_gtf,
            sample_participant_ids = sample_participant_ids,
            vcf_chr_list = vcf_chr_list,
            sample_list = sample_list,
            tpm_threshold = tpm_threshold,
            count_threshold = count_threshold,
            sample_frac_threshold = sample_frac_threshold,
            normalization_method = normalization_method,
            flags = flags
    }

    call compute_PCs {
        input:
            expression_bed = eqtl_prepare_expression.expression_bed,
            prefix = prefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            genotype_covariates = genotype_covariates
    }
}
