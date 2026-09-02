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
    Array[String] cell_types
    Array[String] cell_type_slugs
    Array[File] int_beds
    Array[File] scaled_beds
    Array[File] int_phenotype_pcs
    Array[File] int_phenotype_pcs_all
    Array[File] scaled_phenotype_pcs
    Array[File] scaled_phenotype_pcs_all
    Array[File] int_merged_covariates
    Array[File] scaled_merged_covariates
    Array[File] int_connectivity_outliers
    Array[File] scaled_connectivity_outliers
    String docker_image
    Int cpu = 1
    String memory = "4 GB"
    Int disk_gb = 20
    Int preemptible_attempts = 0
    Int max_retries = 1
  }

  command <<<
    set -euo pipefail
    stage="build_qtl_manifest"
    log="$stage.log"
    status=0
    printf 'stage=%s start_time=%s\n' "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"
    trap 'status=$?; printf "stage=%s error_status=%s time=%s\\n" "$stage" "$status" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$log"; exit "$status"' ERR
    mkdir -p inputs outputs
    cat > inputs/cell_types.txt <<'CELL_TYPES'
~{sep='\n' cell_types}
CELL_TYPES
    cat > inputs/cell_type_slugs.txt <<'CELL_TYPE_SLUGS'
~{sep='\n' cell_type_slugs}
CELL_TYPE_SLUGS
    cat > inputs/int_beds.txt <<'INT_BEDS'
~{sep='\n' int_beds}
INT_BEDS
    cat > inputs/scaled_beds.txt <<'SCALED_BEDS'
~{sep='\n' scaled_beds}
SCALED_BEDS
    cat > inputs/int_phenotype_pcs.txt <<'INT_PHENOTYPE_PCS'
~{sep='\n' int_phenotype_pcs}
INT_PHENOTYPE_PCS
    cat > inputs/int_phenotype_pcs_all.txt <<'INT_PHENOTYPE_PCS_ALL'
~{sep='\n' int_phenotype_pcs_all}
INT_PHENOTYPE_PCS_ALL
    cat > inputs/scaled_phenotype_pcs.txt <<'SCALED_PHENOTYPE_PCS'
~{sep='\n' scaled_phenotype_pcs}
SCALED_PHENOTYPE_PCS
    cat > inputs/scaled_phenotype_pcs_all.txt <<'SCALED_PHENOTYPE_PCS_ALL'
~{sep='\n' scaled_phenotype_pcs_all}
SCALED_PHENOTYPE_PCS_ALL
    cat > inputs/int_merged_covariates.txt <<'INT_MERGED_COVARIATES'
~{sep='\n' int_merged_covariates}
INT_MERGED_COVARIATES
    cat > inputs/scaled_merged_covariates.txt <<'SCALED_MERGED_COVARIATES'
~{sep='\n' scaled_merged_covariates}
SCALED_MERGED_COVARIATES
    cat > inputs/int_connectivity_outliers.txt <<'INT_CONNECTIVITY_OUTLIERS'
~{sep='\n' int_connectivity_outliers}
INT_CONNECTIVITY_OUTLIERS
    cat > inputs/scaled_connectivity_outliers.txt <<'SCALED_CONNECTIVITY_OUTLIERS'
~{sep='\n' scaled_connectivity_outliers}
SCALED_CONNECTIVITY_OUTLIERS
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/build_qtl_manifest.R \
      --cell-types inputs/cell_types.txt \
      --cell-type-slugs inputs/cell_type_slugs.txt \
      --int-beds inputs/int_beds.txt \
      --scaled-beds inputs/scaled_beds.txt \
      --int-pcs inputs/int_phenotype_pcs.txt \
      --int-pcs-all inputs/int_phenotype_pcs_all.txt \
      --scaled-pcs inputs/scaled_phenotype_pcs.txt \
      --scaled-pcs-all inputs/scaled_phenotype_pcs_all.txt \
      --int-covariates inputs/int_merged_covariates.txt \
      --scaled-covariates inputs/scaled_merged_covariates.txt \
      --int-outliers inputs/int_connectivity_outliers.txt \
      --scaled-outliers inputs/scaled_connectivity_outliers.txt \
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
