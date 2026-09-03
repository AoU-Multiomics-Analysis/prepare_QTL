version 1.0

task PrepareHspeBatches {
  input {
    File expression
    Float log2_pseudocount
    File gtf
    File lm22
    Float min_overlap = 0.80
    Float marker_fraction = 0.10
    Int batch_size = 100
    Boolean quantile_normalize = false
    Int random_seed = 20260901
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 2
    Int max_retries = 2
  }

  String quantile_normalize_argument = if quantile_normalize then "--quantile-normalize" else ""

  command <<<
    set -euo pipefail
    stage="prepare_hspe_batches"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/prepare_hspe_batches.R \
      --expression '~{expression}' \
      --gtf '~{gtf}' \
      --lm22 '~{lm22}' \
      --min-overlap '~{min_overlap}' \
      --marker-fraction '~{marker_fraction}' \
      --batch-size '~{batch_size}' \
      --random-seed '~{random_seed}' \
      ~{quantile_normalize_argument} \
      --log2-pseudocount '~{log2_pseudocount}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/hspe_markers.tsv)" \
      "batches,prepared,markers,overlap_report,transformed_lm22" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[File] batches = glob("outputs/batch_*.rds")
    File prepared = "outputs/hspe_prepared.rds"
    File markers = "outputs/hspe_markers.tsv"
    File overlap_report = "outputs/hspe_overlap.tsv"
    File transformed_lm22 = "outputs/hspe_lm22_log.tsv.gz"
    File log = "prepare_hspe_batches.log"
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

task RunHspeBatch {
  input {
    File prepared
    File batch
    String docker_image
    String memory = "4 GB"
    Int disk_gb = 10
    Int preemptible_attempts = 2
    Int max_retries = 2
  }
  command <<<
    set -euo pipefail
    stage="run_hspe_batch"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/run_hspe_batch.R \
      '~{prepared}' '~{batch}' outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=one_batch outputs=hspe_batch_result.rds completion_time=%s\n' \
      "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>
  output {
    File result = "outputs/hspe_batch_result.rds"
    File log = "run_hspe_batch.log"
  }
  runtime {
    docker: docker_image
    cpu: 1
    memory: memory
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: max_retries
  }
}

task MergeHspeBatches {
  input {
    File prepared
    Array[File] batch_results
    File preparation_log
    Array[File] batch_logs
    String docker_image
    Int cpu = 2
    String memory = "16 GB"
    Int disk_gb = 50
    Int preemptible_attempts = 2
    Int max_retries = 2
  }
  command <<<
    set -euo pipefail
    stage="merge_hspe_batches"
    log="run_hspe.log"
    status=0
    tee -a "$log" < '~{preparation_log}'
    while IFS= read -r batch_log; do
      tee -a "$log" < "$batch_log"
    done < '~{write_lines(batch_logs)}'
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/merge_hspe_batches.R \
      '~{prepared}' '~{write_lines(batch_results)}' outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=proportions,metadata,diagnostics completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/hspe_proportions.tsv)" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>
  output {
    File proportions = "outputs/hspe_proportions.tsv"
    File metadata = "outputs/hspe_metadata.json"
    File diagnostics = "outputs/hspe_sample_diagnostics.tsv"
    File log = "run_hspe.log"
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
