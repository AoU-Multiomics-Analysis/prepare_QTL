version 1.0

task FilterCellTypeBed {
  input {
    File cell_type_bed_inventory
    File cell_type_bed
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
  command <<<
    set -euo pipefail
    stage="filter_cell_type_bed"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s status=failed error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    inventory_path='~{sub(cell_type_bed_inventory, "'", "'\"'\"'")}'
    bed_path='~{sub(cell_type_bed, "'", "'\"'\"'")}'
    reference_path='~{if defined(reference_summary) then sub(select_first([reference_summary]), "'", "'\"'\"'") else ""}'
    optional_arguments=()
    if [[ -n "$reference_path" ]]; then
      optional_arguments+=(--reference-summary "$reference_path")
    fi
    residual_value='~{default="" residual_cutoff}'
    if [[ -n "$residual_value" ]]; then
      optional_arguments+=(--residual-cutoff "$residual_value")
    fi
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/downstream/filter_cell_type_beds.R \
      --inventory "$inventory_path" --single-bed "$bed_path" \
      --min-mean-log2-cpm1 '~{min_mean_log2_cpm1}' "${optional_arguments[@]}" \
      --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=cell_types:1 outputs=filtered_bed,reports completion_time=%s\n' \
      "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>
  output {
    File filtered_bed = glob("outputs/beds/*.filtered.bed.gz")[0]
    File filtered_inventory = "outputs/filtered_inventory.tsv"
    File gene_comparison = "outputs/gene_comparison.tsv.gz"
    File filter_metrics = "outputs/filter_metrics.tsv"
    File sample_ids = "outputs/sample_ids.txt"
    File log = "filter_cell_type_bed.log"
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

task MergeFilterReports {
  input {
    File cell_type_bed_inventory
    Array[File] inventories
    Array[File] comparisons
    Array[File] metrics
    Array[File] samples
    Array[File] logs
    Boolean post_residual = false
    String docker_image
    Int cpu = 1
    String memory = "8 GB"
    Int disk_gb = 20
    Int preemptible_attempts = 2
    Int max_retries = 2
  }
  command <<<
    set -euo pipefail
    stage="merge_filter_reports"
    log="$stage.log"
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s status=failed error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    # File arrays are localized before these task-local lists are written.
    cp '~{write_lines(inventories)}' inventories.txt
    cp '~{write_lines(comparisons)}' comparisons.txt
    cp '~{write_lines(metrics)}' metrics.txt
    cp '~{write_lines(samples)}' samples.txt
    cp '~{write_lines(logs)}' logs.txt
    inventory_path='~{sub(cell_type_bed_inventory, "'", "'\"'\"'")}'
    optional_arguments=()
    if [[ '~{post_residual}' == true ]]; then
      optional_arguments+=(--post-residual)
    fi
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/downstream/merge_filter_reports.R \
      --inventory "$inventory_path" --inventories inventories.txt --comparisons comparisons.txt \
      --metrics metrics.txt --samples samples.txt --logs logs.txt \
      "${optional_arguments[@]}" --output-dir outputs 2>&1 | tee -a "$log"
    printf 'stage=%s dimensions=cell_types:%s outputs=filtered_inventory,reports,plots completion_time=%s\n' \
      "$stage" '~{length(inventories)}' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>
  output {
    File filtered_inventory = "outputs/filtered_inventory.tsv"
    File negative_summary = "outputs/negative_summary.tsv.gz"
    File gene_comparison = "outputs/gene_comparison.tsv.gz"
    File filter_metrics = "outputs/filter_metrics.tsv"
    Array[File] plots = glob("outputs/plots/*.pdf")
    File log = "merge_filter_reports.log"
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
