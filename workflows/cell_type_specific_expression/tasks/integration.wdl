version 1.0

task PrepareScatterInputs {
  input {
    File cell_type_bed_inventory
    Array[File] cell_type_beds
    String output_prefix
    String docker_image
    Int cpu = 1
    String memory = "4 GB"
    Int disk_gb = 20
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  File output_prefix_file = write_lines([output_prefix])

  command <<<
    set -euo pipefail
    stage="prepare_scatter_inputs"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p scatter
    cat > scatter/cell_type_bed_paths.txt <<'CELL_TYPE_BED_PATHS'
~{sep='\n' cell_type_beds}
CELL_TYPE_BED_PATHS
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/prepare_scatter_inputs.R \
      --inventory '~{cell_type_bed_inventory}' \
      --bed-paths scatter/cell_type_bed_paths.txt \
      --output-prefix-file '~{output_prefix_file}' \
      --output-dir scatter \
      --log-file scatter/prepare_scatter_inputs.log 2>&1 | tee -a "$log"
    validated_cell_count="$(wc -l < scatter/cell_types.txt)"
    output_paths="scatter/cell_types.txt,scatter/cell_type_slugs.txt,scatter/expression_beds.txt,scatter/output_prefixes.txt"
    printf 'stage=%s dimensions=validated_cell_count:%s outputs=%s output_paths=%s completion_time=%s\n' \
      "$stage" "$validated_cell_count" "cell_types,cell_type_slugs,expression_beds,output_prefixes" \
      "$output_paths" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    Array[String] cell_types = read_lines("scatter/cell_types.txt")
    Array[String] cell_type_slugs = read_lines("scatter/cell_type_slugs.txt")
    Array[File] expression_beds = cell_type_beds
    Array[String] output_prefixes = read_lines("scatter/output_prefixes.txt")
    File log = "prepare_scatter_inputs.log"
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

task BuildQtlManifest {
  input {
    # Path-only metadata: File -> String coercion at this call boundary keeps
    # the original Cromwell output URLs instead of task-localized paths.
    Array[String] cell_types
    Array[String] cell_type_slugs
    Array[String] int_beds
    Array[String] scaled_beds
    Array[String] int_phenotype_pcs
    Array[String] int_phenotype_pcs_all
    Array[String] scaled_phenotype_pcs
    Array[String] scaled_phenotype_pcs_all
    Array[String] int_merged_covariates
    Array[String] scaled_merged_covariates
    Array[String] int_connectivity_outliers
    Array[String] scaled_connectivity_outliers
    # Inventories are localized for slug alignment. BED and report values stay
    # String metadata so that Cromwell preserves their upstream cloud URLs.
    Array[String] source_beds
    File source_bed_inventory
    Array[String] filtered_beds
    File filtered_bed_inventory
    String negative_summary
    String gene_comparison
    String filter_metrics
    String docker_image
    Int cpu = 1
    String memory = "4 GB"
    Int disk_gb = 20
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  File source_inventory_path_file = write_lines([source_bed_inventory])
  File filtered_inventory_path_file = write_lines([filtered_bed_inventory])

  command <<<
    set -euo pipefail
    stage="build_qtl_manifest"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p outputs
    IFS= read -r source_inventory_path < '~{source_inventory_path_file}'
    IFS= read -r filtered_inventory_path < '~{filtered_inventory_path_file}'
    # Task-scoped serialization preserves quoting and keeps file creation off
    # the Terra workflow engine. These strings are upstream output URLs.
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/build_qtl_manifest.R \
      --cell-types '~{write_json(cell_types)}' \
      --cell-type-slugs '~{write_json(cell_type_slugs)}' \
      --int-beds '~{write_json(int_beds)}' \
      --scaled-beds '~{write_json(scaled_beds)}' \
      --int-pcs '~{write_json(int_phenotype_pcs)}' \
      --int-pcs-all '~{write_json(int_phenotype_pcs_all)}' \
      --scaled-pcs '~{write_json(scaled_phenotype_pcs)}' \
      --scaled-pcs-all '~{write_json(scaled_phenotype_pcs_all)}' \
      --int-covariates '~{write_json(int_merged_covariates)}' \
      --scaled-covariates '~{write_json(scaled_merged_covariates)}' \
      --int-outliers '~{write_json(int_connectivity_outliers)}' \
      --scaled-outliers '~{write_json(scaled_connectivity_outliers)}' \
      --source-beds '~{write_json(source_beds)}' \
      --source-inventory "$source_inventory_path" \
      --filtered-beds '~{write_json(filtered_beds)}' \
      --filtered-inventory "$filtered_inventory_path" \
      --negative-summary '~{write_json(negative_summary)}' \
      --gene-comparison '~{write_json(gene_comparison)}' \
      --filter-metrics '~{write_json(filter_metrics)}' \
      --output outputs/cell_type_qtl_manifest.tsv 2>&1 | tee -a "$log"
    validated_cell_count="$(awk 'END { print NR - 1 }' outputs/cell_type_qtl_manifest.tsv)"
    printf 'stage=%s dimensions=validated_cell_count:%s outputs=%s manifest_path=%s completion_time=%s\n' \
      "$stage" "$validated_cell_count" "cell_type_qtl_manifest" \
      "outputs/cell_type_qtl_manifest.tsv" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
  >>>

  output {
    File manifest = "outputs/cell_type_qtl_manifest.tsv"
    File log = "build_qtl_manifest.log"
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
