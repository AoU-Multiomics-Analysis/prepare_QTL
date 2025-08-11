version 1.0


workflow ComputeGenotypePCs {
    input {
        File VCF 
        String OutputPrefix 
        File genotype_rscript
    }
    
    call RComputeGenotypePCs {
        input:
            vcf_file = VCF,
            output_prefix = OutputPrefix,
            genotype_rscript = genotype_rscript
    }
    
    output {
        File GenotypePCs = RComputeGenotypePCs.output_tsv

    }
}


task RComputeGenotypePCs {
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
            memory: "96G"
            cpu: 2
            disks: "local-disk 100 SSD"
        }
    
        output {
            File output_tsv = "~{output_prefix}_genetic_PCs.tsv"
        }



}
