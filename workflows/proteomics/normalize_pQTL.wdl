version 1.0

task NormalizeProteomics {
    input {
        Array[File] OlinkData
        String OutputPrefix
        String ReferencePlate = "000171002612_A1_01-17-2024_12-21-52"

        Int memory
        Int disk_space
        Int num_threads
        String docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
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
        String normalize_pqtl__normalize_proteomics_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:237c02268a4797c7ec72544a8584b16fc82cb5678958d59eb5cf1647e02b0993"
        Array[File] OlinkData
        String OutputPrefix
        String ReferencePlate = "000171002612_A1_01-17-2024_12-21-52"

        Int memory
        Int disk_space
        Int num_threads

    }

    call NormalizeProteomics {
        input:
            OlinkData = OlinkData,
            OutputPrefix = OutputPrefix,
            ReferencePlate = ReferencePlate,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            docker_image = normalize_pqtl__normalize_proteomics_image
    }

    output {
        File MedianNormalizedData = NormalizeProteomics.MedianNormalizedData
        File FilteredProteomicsData = NormalizeProteomics.FilteredProteomicsData
    }
}
