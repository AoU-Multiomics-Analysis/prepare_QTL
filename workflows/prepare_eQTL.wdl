version 1.0

import "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/prepare_QTL/refs/heads/main/workflows/calculate_phenotypePCs.wdl" as calculate_phenotypePCs


task eqtl_prepare_expression {
    input {
        File tpm_gct
        File counts_gct
        File annotation_gtf
        File sample_participant_ids
        File vcf_chr_list
        File sample_list
        String OutputPrefix 
        

        Int memory 
        Int disk_space 
        Int num_threads 

        Float? tpm_threshold
        Int? count_threshold
        Float? sample_frac_threshold
        String? normalization_method
        String? flags  # --convert_tpm, --legacy_mode
    }
    command {
        set -euo pipefail
        /src/eqtl_prepare_expression.py ${tpm_gct} ${counts_gct} \
        ${annotation_gtf} ${sample_participant_ids} ${vcf_chr_list} ${OutputPrefix} \
        ${"--tpm_threshold " + tpm_threshold} \
        ${"--sample_ids " + sample_list} \
        ${"--count_threshold " + count_threshold} \
        ${"--sample_frac_threshold " + sample_frac_threshold} \
        ${"--normalization_method " + normalization_method} \
        ${flags}
    }

    runtime {
        docker: "quay.io/jonnguye/modified_gtex_eqtl:1.1"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File ExpressionBed="${OutputPrefix}.expression.bed.gz"
        File ExpressionBedIndex="${OutputPrefix}.expression.bed.gz.tbi"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow eQTLPrepareData {
        String OutputPrefix 
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
            OutputPrefix = OutputPrefix,
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

    call calculate_phenotypePCs.ComputePCs {
        input:
            BedFile = eqtl_prepare_expression.ExpressionBed,
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
}
