version 1.0
import  "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/prepare_QTL/refs/heads/main/workflows/calculate_phenotypePCs.wdl" as ComputePCs
import  "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/prepare_QTL/refs/heads/main/workflows/MergeCovariates.wdl" as CovariateMerge




task PrepareProteomicData {
    input {
        File AnnotationGTF
        File SampleList 
        File ProteomicData 
        String OutputPrefix 
        Boolean RankNormalize = true
        
        Int memory 
        Int disk_space 
        Int num_threads
    }
    command {
        Rscript /tmp/PrepareProteomics.R \
            --ProteomicData ${ProteomicData} \
            --AnnotationGTF ${AnnotationGTF} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix} \
            --RankNormalize ${RankNormalize}
        }

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    output {
        File ProteomicBed="${OutputPrefix}.protein.bed.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow pQTLPrepareData {
    input {
        Int memory 
        Int disk_space 
        Int num_threads 
        File AnnotationGTF 
        File SampleList 
        File ProteomicData
        String OutputPrefix 
        Boolean RankNormalize = true
        File? AdditionalCovariates
    } 
    call PrepareProteomicData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix,
            ProteomicData = ProteomicData,
            RankNormalize = RankNormalize

    }

    call ComputePCs.PhenotypePCs {
        input:
            BedFile = PrepareProteomicData.ProteomicBed,
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
    if (defined(AdditionalCovariates)) {
        call CovariateMerge.MergeCovariates as MergeAdditionalCovariates {
            input:
                GenotypePCs = select_first([AdditionalCovariates]),
                MolecularPCs = PhenotypePCs.OutPhenotypePCs,
                OutputPrefix = OutputPrefix
        }
    }
    output {
        File BedFile = PrepareProteomicData.ProteomicBed 
        File PhenotypePCsOut = PhenotypePCs.OutPhenotypePCs 
        File? QtlCovariates = MergeAdditionalCovariates.QtlCovariates
    }
}
