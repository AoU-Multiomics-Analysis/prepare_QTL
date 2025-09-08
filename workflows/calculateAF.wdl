version 1.0 

workflow CaclulateAF {
    input {
        File pvar
        File psam
        File pgen
        String prefix
        Int memory 
        Int disk_space 
        Int num_threads 
        }
    call plink2AF {
        input:
            pvar = pvar, 
            psam = psam,
            pgen = pgen,
            prefix = prefix,
            memory = memory, 
            disk_space = disk_space, 
            num_threads = num_threads
        }
    output {
        File PlinkAF = plink2AF.PlinkAF
        }
    }

    task plink2AF {
        input {
            File pvar 
            File pgen 
            File psam
            String prefix
            Int memory 
            Int disk_space 
            Int num_threads
        }

        command <<< 
        plink2 --pfile ~{prefix} --freq --out ~{prefix}         
        >>> 
        
        runtime {
            docker: "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
            memory: "~{memory}GB"
            disks: "local-disk ~{disk_space} HDD"
            cpu: "~{num_threads}"

        }

        output {
            File PlinkAF = "~{prefix}.afreq" 
        }
    }
