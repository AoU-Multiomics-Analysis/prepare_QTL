version 1.0

# Annotation of cohort-filtered methylation sites.

task AnnotateMethylationSites {
    input {
        File PassingSiteMetadata
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations
        String OutputPrefix
        Int PromoterWindow
        Int MemoryGB
        Int DiskGB
        String docker_image
    }

    command <<<
        Rscript /tmp/AnnotateMethylationSites.R \
            --PassingSiteMetadata "~{PassingSiteMetadata}" \
            --AnnotationGTF "~{AnnotationGTF}" \
            --CCREAnnotations "~{CCREAnnotations}" \
            --CpGIslandAnnotations "~{CpGIslandAnnotations}" \
            --OutputPrefix "~{OutputPrefix}" \
            --PromoterWindow ~{PromoterWindow}
    >>>

    runtime {
        docker: docker_image
        memory: "~{MemoryGB}G"
        disks: "local-disk ~{DiskGB} HDD"
        cpu: 1
    }

    output {
        File PassingSiteAnnotations = "~{OutputPrefix}.methylation.passing_site_annotations.tsv.gz"
    }
}
workflow AnnotateMethylationCohortSites {
    input {
        File PassingSiteMetadata
        File AnnotationGTF
        File CCREAnnotations
        File CpGIslandAnnotations
        String OutputPrefix
        Int PromoterWindow
        Int MemoryGB
        Int DiskGB
        String methylation_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
    }

    call AnnotateMethylationSites {
        input:
            PassingSiteMetadata = PassingSiteMetadata,
            AnnotationGTF = AnnotationGTF,
            CCREAnnotations = CCREAnnotations,
            CpGIslandAnnotations = CpGIslandAnnotations,
            OutputPrefix = OutputPrefix,
            PromoterWindow = PromoterWindow,
            MemoryGB = MemoryGB,
            DiskGB = DiskGB,
            docker_image = methylation_docker_image
    }

    output {
        File PassingSiteAnnotations = AnnotateMethylationSites.PassingSiteAnnotations
    }
}
