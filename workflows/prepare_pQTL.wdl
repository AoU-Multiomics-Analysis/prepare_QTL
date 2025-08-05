version 1.0
import "https://raw.githubusercontent.com/AoU-Multiomics-Analysis/prepare_QTL/refs/heads/main/workflows/calculate_phenotypePCs.wdl" as ComputePCs 




task PrepareProteomicData {
    input {
        File AnnotationGTF
        File SampleList 
        File ProteomicData 
        String OutputPrefix 
        
        Int memory 
        Int disk_space 
        Int num_threads
    }
    command {
        Rscript /tmp/PrepareProteomics.R \
            --ProteomicData ${ProteomicData} \
            --AnnotationGTF ${AnnotationGTF} \
            --SampleList ${SampleList} \
            --OutputPrefix ${OutputPrefix}
        }

    runtime {
        docker: "ghcr.io/AoU-Multiomics-Analysis/prepare_QTL:main"
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
        String OutputPrefix 
    } 
    call PrepareProteomicData {
        input:
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads,
            AnnotationGTF = AnnotationGTF,
            SampleList = SampleList,
            OutputPrefix = OutputPrefix
    }

    call ComputePCs.ComputePCs {
        input:
            BedFile = PrepareProteomicData.ProteomicBed,
            OutputPrefix = OutputPrefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
    output {
        File BedFile = PrepareProteomicData.ProteomicBed 
        File PhenotypePCs = ComputePCs.OutPhenotypePCs 
    }
}
