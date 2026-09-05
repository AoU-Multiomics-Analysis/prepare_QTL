version 1.0

task NormalizeProteomics {
    input {
        Array[File] OlinkData
        String OutputPrefix
        String ReferencePlate = "000171002612_A1_01-17-2024_12-21-52"

        Int memory
        Int disk_space
        Int num_threads
        String docker_image
    }

    command <<<
        mkdir olink_data
        for olink_file in ~{sep=' ' OlinkData}; do
            ln -s "${olink_file}" olink_data/
        done

        Rscript /tmp/NormalizeProteomics.R \
            --OlinkDataDir olink_data \
            --OutputPrefix ~{OutputPrefix} \
            --OutputDir . \
            --ReferencePlate ~{ReferencePlate}
    >>>

    runtime {
        docker: docker_image
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File MedianNormalizedData = "~{OutputPrefix}_median_normalized.tsv.gz"
        File FilteredProteomicsData = "~{OutputPrefix}_npx_values.tsv.gz"
    }
}

workflow NormalizeProteomicsData {
    input {
        Array[File] OlinkData
        String OutputPrefix
        String ReferencePlate = "000171002612_A1_01-17-2024_12-21-52"

        Int memory
        Int disk_space
        Int num_threads
        String proteomics_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845"
    }

    call NormalizeProteomics {
        input:
            OlinkData = OlinkData,
            OutputPrefix = OutputPrefix,
            ReferencePlate = ReferencePlate,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            docker_image = proteomics_docker_image
    }

    output {
        File MedianNormalizedData = NormalizeProteomics.MedianNormalizedData
        File FilteredProteomicsData = NormalizeProteomics.FilteredProteomicsData
    }
}
