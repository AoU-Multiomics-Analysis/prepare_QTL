version 1.0

task PrepareHaemopedia {
  input {
    File counts
    String docker_image
    Int cpu = 1
    String memory = "8 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  File counts_path_file = write_lines([counts])

  command <<<
    set -euo pipefail
    stage="prepare_haemopedia"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s status=failed error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p outputs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    IFS= read -r counts_path < '~{counts_path_file}'
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/prepare_haemopedia.R \
      "$counts_path" outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=reference_prepared outputs=%s completion_time=%s\n' "$stage" \
      'reference_summary.tsv.gz,reference_samples.tsv,reference_metadata.json' \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File summary = "outputs/reference_summary.tsv.gz"
    File samples = "outputs/reference_samples.tsv"
    File metadata = "outputs/reference_metadata.json"
    File log = "prepare_haemopedia.log"
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

task FilterCellTypeBeds {
  input {
    File cell_type_bed_inventory
    Array[File] cell_type_beds
    File? reference_summary
    Float min_mean_log2_cpm1 = 0.01
    Float? residual_cutoff
    String docker_image
    Int cpu = 1
    String memory = "8 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  File config_json = write_json(object {
    inventory: cell_type_bed_inventory,
    bed_paths: cell_type_beds,
    reference_summary: reference_summary,
    min_mean_log2_cpm1: min_mean_log2_cpm1,
    residual_cutoff: residual_cutoff
  })

  command <<<
    set -euo pipefail
    stage="filter_cell_type_beds"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s status=failed error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p outputs
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/filter_cell_type_beds.R \
      '~{config_json}' outputs 2>&1 | tee -a "$log"
    retained_count="$(awk 'END { print NR - 1 }' outputs/filtered_inventory.tsv)"
    printf 'stage=%s dimensions=cell_types:%s outputs=%s completion_time=%s\n' "$stage" \
      "$retained_count" 'filtered_beds,filtered_inventory,negative_summary,gene_comparison,filter_metrics,plots' \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] filtered_beds = read_lines("outputs/filtered_beds.txt")
    File filtered_inventory = "outputs/filtered_inventory.tsv"
    File negative_summary = "outputs/negative_summary.tsv.gz"
    File gene_comparison = "outputs/gene_comparison.tsv.gz"
    File filter_metrics = "outputs/filter_metrics.tsv"
    Array[File] plots = glob("outputs/plots/*.pdf")
    File log = "filter_cell_type_beds.log"
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
