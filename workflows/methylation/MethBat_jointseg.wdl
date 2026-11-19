version 1.0

workflow MethBat_jointseg {
    meta {
        description: "Cohort joint segmentation with MethBat v0.17 joint-segment. Averages pb-CpG-tools methylation across samples and writes Methylated, Unmethylated, and AlleleSpecificMethylation regions."
    }

    input {
        # TSV, optional header. Columns match the v0.17 collection file:
        #   identifier, filename, labels
        # filename is a pb-CpG-tools output prefix (the stem of prefix.combined.bed.gz,
        # prefix.hap1.bed.gz, and prefix.hap2.bed.gz), not a single BED path.
        # labels may be blank or ';'-delimited and are copied through unchanged.
        File manifest_file
        String prefix

        Int threads = 8
        # joint-segment loads every dataset's CpGs into memory (~1 GB RAM per WGS sample).
        Int mem_gb = 64
        # Disk must hold the localized combined/hap1/hap2 BEDs plus the segment outputs.
        Int disk_gb = 1000
        String? gcs_billing_project

        Boolean condense_bed_labels = false
        Boolean enable_nodata_segments = false

        RuntimeAttr? runtime_attr_override
    }

    call MethBat_joint_segmentation {
        input:
            manifest_file = manifest_file,
            prefix = prefix,
            threads = threads,
            mem_gb = mem_gb,
            disk_gb = disk_gb,
            gcs_billing_project = gcs_billing_project,
            condense_bed_labels = condense_bed_labels,
            enable_nodata_segments = enable_nodata_segments,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File meth_regions_bed = MethBat_joint_segmentation.meth_regions_bed
        File combined_methyl_bedgraph = MethBat_joint_segmentation.combined_methyl_bedgraph
        File asm_bedgraph = MethBat_joint_segmentation.asm_bedgraph
    }
}


task MethBat_joint_segmentation {
    meta {
        description: "Build a v0.17 methbat input-collection from pb-CpG-tools prefixes and run methbat joint-segment."
    }

    input {
        File manifest_file
        String prefix
        Int threads
        Int mem_gb
        Int disk_gb
        String? gcs_billing_project
        Boolean condense_bed_labels
        Boolean enable_nodata_segments
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p output input_cpg
        cohort_tsv="output/~{prefix}.cohort.tsv"
        billing="~{if defined(gcs_billing_project) then gcs_billing_project else ""}"

        # Drop a header row and blank lines. labels (column 3) may be empty.
        tr -d '\r' < "~{manifest_file}" | \
            awk -F '\t' 'NR==1 && $1=="identifier" {next} $0 ~ /^[[:space:]]*$/ {next} {print}' \
            > manifest.body.tsv

        cut -f 1 manifest.body.tsv > sampleid.txt
        cut -f 2 manifest.body.tsv > prefixlist.txt
        awk -F '\t' '{print $3}' manifest.body.tsv > labels.txt

        n_prefixes=$(wc -l < prefixlist.txt | tr -d ' ')
        n_samples=$(wc -l < sampleid.txt | tr -d ' ')
        n_labels=$(wc -l < labels.txt | tr -d ' ')
        if [[ "${n_prefixes}" -eq 0 ]]; then
            echo "manifest has no sample rows" >&2
            exit 1
        fi
        if [[ "${n_prefixes}" -ne "${n_samples}" || "${n_prefixes}" -ne "${n_labels}" ]]; then
            echo "identifier, filename, and labels columns must have the same number of rows" >&2
            exit 1
        fi

        dups=$(cut -f 1 sampleid.txt | sort | uniq -d | paste -sd, -)
        if [[ -n "${dups}" ]]; then
            echo "duplicate identifiers: ${dups}" >&2
            exit 1
        fi

        printf "identifier\tfilename\tlabels\n" > "${cohort_tsv}"

        idx=0
        while IFS=$'\t' read -r identifier cpg_prefix label || [[ -n "${identifier}" ]]; do
            if [[ -z "${identifier}" || -z "${cpg_prefix}" ]]; then
                echo "each manifest row needs an identifier and a pb-CpG-tools prefix" >&2
                exit 1
            fi

            idx=$((idx + 1))
            printf -v n "%05d" "${idx}"
            dest="input_cpg/${n}"
            mkdir -p "${dest}"
            # filename is a prefix: copy prefix.combined/hap1/hap2.bed.gz and sibling indexes.
            if [[ "${cpg_prefix}" == gs://* ]]; then
                if [[ -n "${billing}" ]]; then
                    gsutil -u "${billing}" cp "${cpg_prefix}*" "${dest}/"
                else
                    gsutil cp "${cpg_prefix}*" "${dest}/"
                fi
            else
                matches=( "${cpg_prefix}"* )
                if [[ ${#matches[@]} -eq 0 || ! -e "${matches[0]}" ]]; then
                    echo "no pb-CpG-tools files at prefix ${cpg_prefix}" >&2
                    exit 1
                fi
                cp "${matches[@]}" "${dest}/"
            fi

            local_prefix="${dest}/$(basename "${cpg_prefix}")"
            for suffix in combined hap1 hap2; do
                if [[ ! -f "${local_prefix}.${suffix}.bed.gz" ]]; then
                    echo "missing ${local_prefix}.${suffix}.bed.gz (filename must be a pb-CpG-tools prefix)" >&2
                    exit 1
                fi
            done
            printf "%s\t%s\t%s\n" "${identifier}" "${local_prefix}" "${label}" >> "${cohort_tsv}"
        done < <(paste sampleid.txt prefixlist.txt labels.txt)

        extra_args=()
        if [[ "~{condense_bed_labels}" == "true" ]]; then
            extra_args+=(--condense-bed-labels)
        fi
        if [[ "~{enable_nodata_segments}" == "true" ]]; then
            extra_args+=(--enable-nodata-segments)
        fi

        if [[ ${#extra_args[@]} -eq 0 ]]; then
            methbat joint-segment \
                -t ~{threads} \
                --input-collection "${cohort_tsv}" \
                --output-prefix "output/~{prefix}"
        else
            methbat joint-segment \
                -t ~{threads} \
                --input-collection "${cohort_tsv}" \
                --output-prefix "output/~{prefix}" \
                "${extra_args[@]}"
        fi
    >>>

    output {
        File meth_regions_bed = "output/~{prefix}.meth_regions.bed"
        File combined_methyl_bedgraph = "output/~{prefix}.combined_methyl.bedgraph"
        File asm_bedgraph = "output/~{prefix}.asm.bedgraph"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             mem_gb,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}