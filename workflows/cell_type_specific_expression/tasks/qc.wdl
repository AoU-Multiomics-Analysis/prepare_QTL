version 1.0

task BuildManifest {
  input {
    File cell_type_bed_inventory
    File export_qc_summary
    File export_qc_plots
    File model
    File model_log
    File original_proportions
    File combined_proportions
    File tca_weights
    File filter_report
    File? hspe_metadata
    String proportion_mode
    Float log2_pseudocount
    Float min_lm22_overlap
    Float hspe_marker_fraction
    String hspe_marker_method
    Boolean hspe_quantile_normalize
    Float group_mean_threshold
    Float zero_floor
    Int tca_max_iters
    Boolean tca_parallel
    Array[String] gene_type
    Int random_seed
    String scale
    String tca_version = "1.2.1"
    String container_image
    String docker_image
    Int cpu = 4
    String memory = "32 GB"
    Int disk_gb = 100
    Int preemptible_attempts = 1
    Int max_retries = 2
  }

  String hspe_metadata_path = if defined(hspe_metadata) then select_first([hspe_metadata]) else ""

  command <<<
    set -euo pipefail
    stage="build_manifest"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    public_inventory=outputs/output_inventory.tsv
    hspe_metadata_path="~{hspe_metadata_path}"
    hspe_metadata_arguments=()
    if [[ -n "$hspe_metadata_path" ]]; then
      hspe_metadata_arguments=(--hspe-metadata "$hspe_metadata_path")
    fi
    mkdir -p outputs
    cp '~{cell_type_bed_inventory}' "$public_inventory"
    printf '%s\n' \
      '~{export_qc_summary}' \
      '~{export_qc_plots}' \
      '~{model}' \
      '~{model_log}' \
      '~{original_proportions}' \
      '~{combined_proportions}' \
      '~{tca_weights}' \
      '~{filter_report}' > supporting_inputs.txt
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/build_deconvolution_manifest.R \
      --outputs "$public_inventory" \
      --export-qc-summary '~{export_qc_summary}' \
      --original-proportions '~{original_proportions}' \
      --combined-proportions '~{combined_proportions}' \
      --tca-weights '~{tca_weights}' \
      --filter-report '~{filter_report}' \
      --model '~{model}' \
      --model-log '~{model_log}' \
      "${hspe_metadata_arguments[@]}" \
      --tca-version '~{tca_version}' \
      --proportion-mode '~{proportion_mode}' \
      --log2-pseudocount '~{log2_pseudocount}' \
      --min-lm22-overlap '~{min_lm22_overlap}' \
      --hspe-marker-fraction '~{hspe_marker_fraction}' \
      --hspe-marker-method '~{hspe_marker_method}' \
      --hspe-quantile-normalize '~{hspe_quantile_normalize}' \
      --group-mean-threshold '~{group_mean_threshold}' \
      --zero-floor '~{zero_floor}' \
      --tca-max-iters '~{tca_max_iters}' \
      --tca-parallel '~{tca_parallel}' \
      --gene-types '~{sep="," gene_type}' \
      --random-seed '~{random_seed}' \
      --scale '~{scale}' \
      --effective-parameters-output outputs/effective_parameters.json \
      --container-image '~{container_image}' \
      --output outputs/output_manifest.json \
      --qc-output outputs/qc_summary.tsv \
      --log-file outputs/output_manifest.log 2>&1 | tee -a "$log"
    bed_count="$(awk 'END { print NR - 1 }' "$public_inventory")"
    printf 'stage=%s dimensions=beds:%s outputs=%s completion_time=%s\n' \
      "$stage" "$bed_count" \
      "output_manifest,qc_summary,qc_plots,provenance,effective_parameters" \
      "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File output_manifest = "outputs/output_manifest.json"
    File qc_summary = "outputs/qc_summary.tsv"
    File qc_plots = export_qc_plots
    File provenance = "outputs/output_inventory.tsv"
    File effective_parameters_file = "outputs/effective_parameters.json"
    File log = "build_manifest.log"
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
