version 1.0

task FilterExpressionGenes {
  input {
    File input_expression
    File gtf
    Array[String] gene_type
    Float log2_pseudocount
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 1
    Int max_retries = 2
  }

  command <<<
    set -euo pipefail
    stage="filter_expression_genes"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s gene_types=%s\n' \
      "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" '~{sep="," gene_type}' | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p outputs
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/filter_expression_genes.R \
      --expression '~{input_expression}' \
      --gtf '~{gtf}' \
      --gene-types '~{sep="," gene_type}' \
      --log2-pseudocount '~{log2_pseudocount}' \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s outputs=filtered_expression,gene_type_filter_report completion_time=%s\n' \
      "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File expression = "outputs/filtered_expression.bed.gz"
    File report = "outputs/gene_type_filter_report.tsv"
    File log = "filter_expression_genes.log"
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
