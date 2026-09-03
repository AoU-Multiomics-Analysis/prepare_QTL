version 1.0

task SummarizeCellTypeBeds {
  input {
    Array[File] cell_type_beds
    File cell_type_bed_inventory
    String docker_image
    Int cpu = 1
    String memory = "8 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  command <<<
    set -euo pipefail
    stage="summarize_cell_type_beds"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p outputs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/summarize_cell_type_beds.R \
      '~{cell_type_bed_inventory}' \
      '~{write_lines(cell_type_beds)}' \
      outputs/cell_type_gene_summary.tsv.gz 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=cell_types:%s outputs=%s completion_time=%s\n' \
      "$stage" '~{length(cell_type_beds)}' 'cell_type_gene_summary.tsv.gz' \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File summary = "outputs/cell_type_gene_summary.tsv.gz"
    File log = "summarize_cell_type_beds.log"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: memory
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}
