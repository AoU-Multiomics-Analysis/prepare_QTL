version 1.0

task FitTca {
  input {
    File expression
    Float log2_pseudocount
    File tca_weights
    File? covariates
    Int max_iters = 10
    Int random_seed = 20260901
    Boolean parallel = false
    String docker_image
    Int cpu = 16
    String memory = "256 GB"
    Int disk_gb = 750
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_path = if defined(covariates) then select_first([covariates]) else ""
  String parallel_argument = if parallel then "--parallel" else ""

  command <<<
    set -euo pipefail
    stage="fit_tca"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    covariates_path="~{covariates_path}"
    covariates_arguments=()
    if [[ -n "$covariates_path" ]]; then
      covariates_arguments=(--covariates "$covariates_path")
    fi
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/fit_tca.R \
      --expression '~{expression}' \
      --weights '~{tca_weights}' \
      "${covariates_arguments[@]}" \
      --num-cores '~{cpu}' \
      ~{parallel_argument} \
      --max-iters '~{max_iters}' \
      --random-seed '~{random_seed}' \
      --log2-pseudocount '~{log2_pseudocount}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s outputs=%s completion_time=%s\n' \
      "$stage" "model,model_log,excluded_genes" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File model = "outputs/tca_model.rds"
    File model_log = "outputs/tca_model.log"
    File excluded_genes = "outputs/tca_excluded_genes.tsv"
    File log = "fit_tca.log"
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

task ExportTcaBeds {
  input {
    File expression
    Float log2_pseudocount
    File model
    File tca_weights
    File? covariates
    Boolean parallel = false
    String docker_image
    Int cpu = 8
    String memory = "256 GB"
    Int disk_gb = 500
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  String covariates_path = if defined(covariates) then select_first([covariates]) else ""
  String parallel_argument = if parallel then "--parallel" else ""

  command <<<
    set -euo pipefail
    stage="export_tca_beds"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    covariates_path="~{covariates_path}"
    covariates_arguments=()
    if [[ -n "$covariates_path" ]]; then
      covariates_arguments=(--covariates "$covariates_path")
    fi
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/export_tca_beds.R \
      --expression '~{expression}' \
      --model '~{model}' \
      --weights '~{tca_weights}' \
      "${covariates_arguments[@]}" \
      --num-cores '~{cpu}' \
      ~{parallel_argument} \
      --log2-pseudocount '~{log2_pseudocount}' \
      --output-dir outputs \
      --log-file outputs/export_tca_beds.log 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/cell_type_bed_inventory.tsv)" \
      "cell_type_beds,ordered_bed_paths,inventory,reconstruction,qc_summary,qc_plots" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] cell_type_beds = read_lines("outputs/cell_type_bed_paths.txt")
    File cell_type_bed_paths = "outputs/cell_type_bed_paths.txt"
    File cell_type_bed_inventory = "outputs/cell_type_bed_inventory.tsv"
    File reconstruction_by_sample = "outputs/reconstruction_by_sample.tsv"
    File qc_summary = "outputs/qc_summary.tsv"
    File qc_plots = "outputs/qc_plots.pdf"
    File export_detail_log = "outputs/export_tca_beds.log"
    File log = "export_tca_beds.log"
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
