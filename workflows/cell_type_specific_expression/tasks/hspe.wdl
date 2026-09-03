version 1.0

task RunHspe {
  input {
    File expression
    Float log2_pseudocount
    File gtf
    File lm22
    Float min_overlap = 0.80
    Float marker_fraction = 0.10
    String marker_method = "ratio"
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
    stage="run_hspe"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/run_hspe.R \
      --expression '~{expression}' \
      --gtf '~{gtf}' \
      --lm22 '~{lm22}' \
      --min-overlap '~{min_overlap}' \
      --marker-fraction '~{marker_fraction}' \
      --marker-method '~{marker_method}' \
      --random-seed '~{random_seed}' \
      ~{quantile_normalize_argument} \
      --log2-pseudocount '~{log2_pseudocount}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=%s outputs=%s completion_time=%s\n' \
      "$stage" "$(wc -l < outputs/hspe_proportions.tsv)" \
      "proportions,markers,metadata,overlap_report,transformed_lm22" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File proportions = "outputs/hspe_proportions.tsv"
    File markers = "outputs/hspe_markers.tsv"
    File metadata = "outputs/hspe_metadata.json"
    File overlap_report = "outputs/hspe_overlap.tsv"
    File transformed_lm22 = "outputs/hspe_lm22_log.tsv.gz"
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
