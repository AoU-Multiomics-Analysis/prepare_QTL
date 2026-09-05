version 1.0

import "tasks/hspe.wdl" as hspe_tasks
import "tasks/expression.wdl" as expression_tasks
import "tasks/proportions.wdl" as proportion_tasks
import "tasks/tca.wdl" as tca_tasks
import "tasks/qc.wdl" as qc_tasks
import "tasks/gene_summary.wdl" as summary_tasks
import "tasks/reference_filter.wdl" as reference_filter_tasks

workflow CellTypeDeconvolution {
  input {
    File expression
    File gtf
    File lm22
    File? precomputed_tca_model
    File? precomputed_proportions
    File? covariates
    String estimation_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:9f7af7c16fa3dc7a0b82c042a40145fa26afce4a96547791e0b29a9e8de4d754"
    String fit_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:9f7af7c16fa3dc7a0b82c042a40145fa26afce4a96547791e0b29a9e8de4d754"
    String export_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:9f7af7c16fa3dc7a0b82c042a40145fa26afce4a96547791e0b29a9e8de4d754"
    String downstream_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:9f7af7c16fa3dc7a0b82c042a40145fa26afce4a96547791e0b29a9e8de4d754"
    Int preemptible_attempts = 2
    Int max_retries = 2
    Float min_lm22_overlap = 0.80
    Float hspe_marker_fraction = 0.10
    Boolean hspe_quantile_normalize = false
    Int hspe_batch_size = 100
    String hspe_batch_memory = "4 GB"
    Int hspe_batch_disk_gb = 10
    Float group_mean_threshold = 0.0001
    Float zero_floor = 0.000001
    Int tca_max_iters = 10
    Boolean tca_parallel = false
    Int random_seed = 20260901
    Float log2_pseudocount = 0.0
    Array[String] gene_type = ["protein_coding", "lncRNA"]
    File? haemopedia_counts
    Float reference_min_mean_log2_cpm1 = 0.01
    Float? reference_residual_cutoff

    Int hspe_cpu = 4
    String hspe_memory = "32 GB"
    Int hspe_disk_gb = 100
    Int proportions_cpu = 2
    String proportions_memory = "16 GB"
    Int proportions_disk_gb = 50
    Int fit_cpu = 16
    String fit_memory = "256 GB"
    Int fit_disk_gb = 750
    Int export_cpu = 8
    String export_memory = "256 GB"
    Int export_disk_gb = 500
    String gene_summary_memory = "8 GB"
    Int manifest_cpu = 4
    String manifest_memory = "32 GB"
    Int manifest_disk_gb = 100
  }

  String tca_version = "1.2.1"
  if (!defined(precomputed_tca_model)) {
    File? fit_covariates = covariates
    call expression_tasks.FilterExpressionGenes {
      input:
        input_expression = expression,
        gtf = gtf,
        gene_type = gene_type,
        log2_pseudocount = log2_pseudocount,
        docker_image = estimation_docker_image,
        cpu = hspe_cpu,
        memory = hspe_memory,
        disk_gb = hspe_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }

    call proportion_tasks.ValidateProportionMode {
      input:
        precomputed_proportions = precomputed_proportions,
        docker_image = estimation_docker_image,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }

  }

  String proportion_mode = if defined(precomputed_tca_model) then "precomputed_model"
    else select_first([ValidateProportionMode.selected_mode])
  String hspe_marker_method = "ratio"

  if (!defined(precomputed_tca_model) && !defined(precomputed_proportions)) {
    call hspe_tasks.PrepareHspeBatches {
      input:
        expression = select_first([FilterExpressionGenes.expression]),
        log2_pseudocount = log2_pseudocount,
        gtf = gtf,
        lm22 = lm22,
        min_overlap = min_lm22_overlap,
        marker_fraction = hspe_marker_fraction,
        batch_size = hspe_batch_size,
        quantile_normalize = hspe_quantile_normalize,
        random_seed = random_seed,
        docker_image = estimation_docker_image,
        cpu = hspe_cpu,
        memory = hspe_memory,
        disk_gb = hspe_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
    scatter (batch in PrepareHspeBatches.batches) {
      call hspe_tasks.RunHspeBatch {
        input:
          prepared = PrepareHspeBatches.prepared,
          batch = batch,
          docker_image = estimation_docker_image,
          memory = hspe_batch_memory,
          disk_gb = hspe_batch_disk_gb,
          preemptible_attempts = preemptible_attempts,
          max_retries = max_retries
      }
    }
    call hspe_tasks.MergeHspeBatches {
      input:
        prepared = PrepareHspeBatches.prepared,
        batch_results = RunHspeBatch.result,
        preparation_log = PrepareHspeBatches.log,
        batch_logs = RunHspeBatch.log,
        docker_image = estimation_docker_image,
        cpu = proportions_cpu,
        memory = proportions_memory,
        disk_gb = proportions_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  if (!defined(precomputed_tca_model)) {
    File proportions_for_processing = if !defined(precomputed_proportions)
      then select_first([MergeHspeBatches.proportions])
      else select_first([precomputed_proportions])

    call proportion_tasks.ProcessProportions {
      input:
        proportions = proportions_for_processing,
        mean_threshold = group_mean_threshold,
        zero_floor = zero_floor,
        docker_image = estimation_docker_image,
        cpu = proportions_cpu,
        memory = proportions_memory,
        disk_gb = proportions_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }

    call tca_tasks.FitTca {
      input:
        expression = select_first([FilterExpressionGenes.expression]),
        tca_weights = ProcessProportions.tca_weights,
        covariates = covariates,
        max_iters = tca_max_iters,
        random_seed = random_seed,
        parallel = tca_parallel,
        docker_image = fit_docker_image,
        cpu = fit_cpu,
        memory = fit_memory,
        disk_gb = fit_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }

  }

  File export_expression = if defined(precomputed_tca_model) then expression
    else select_first([FilterExpressionGenes.expression])

  call tca_tasks.CleanTcaModel {
    input:
      unfiltered_model = select_first([precomputed_tca_model, FitTca.model]),
      reuse_model = defined(precomputed_tca_model),
      docker_image = fit_docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call tca_tasks.ExportTcaBeds {
    input:
      expression = export_expression,
      model = CleanTcaModel.model,
      tca_weights = CleanTcaModel.tca_weights,
      covariates = fit_covariates,
      reuse_model = defined(precomputed_tca_model),
      parallel = tca_parallel,
      docker_image = export_docker_image,
      cpu = export_cpu,
      memory = export_memory,
      disk_gb = export_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call summary_tasks.SummarizeCellTypeBeds {
    input:
      cell_type_beds = ExportTcaBeds.cell_type_beds,
      cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory,
      docker_image = downstream_docker_image,
      memory = gene_summary_memory,
      disk_gb = export_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  if (defined(haemopedia_counts)) {
    call reference_filter_tasks.PrepareHaemopedia {
      input:
        counts = select_first([haemopedia_counts]),
        docker_image = downstream_docker_image,
        memory = gene_summary_memory,
        disk_gb = export_disk_gb,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries
    }
  }

  call reference_filter_tasks.FilterCellTypeBeds {
    input:
      cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory,
      cell_type_beds = ExportTcaBeds.cell_type_beds,
      reference_summary = PrepareHaemopedia.summary,
      min_mean_log2_cpm1 = reference_min_mean_log2_cpm1,
      residual_cutoff = reference_residual_cutoff,
      docker_image = downstream_docker_image,
      memory = gene_summary_memory,
      disk_gb = export_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  call qc_tasks.BuildManifest {
    input:
      cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory,
      export_qc_summary = ExportTcaBeds.qc_summary,
      export_qc_plots = ExportTcaBeds.qc_plots,
      model = CleanTcaModel.model,
      model_log = FitTca.model_log,
      original_proportions = ProcessProportions.original,
      combined_proportions = ProcessProportions.combined,
      tca_weights = CleanTcaModel.tca_weights,
      filter_report = ProcessProportions.filter_report,
      hspe_metadata = MergeHspeBatches.metadata,
      proportion_mode = proportion_mode,
      log2_pseudocount = log2_pseudocount,
      min_lm22_overlap = min_lm22_overlap,
      hspe_marker_fraction = hspe_marker_fraction,
      hspe_marker_method = hspe_marker_method,
      hspe_quantile_normalize = hspe_quantile_normalize,
      group_mean_threshold = group_mean_threshold,
      zero_floor = zero_floor,
      tca_max_iters = tca_max_iters,
      tca_parallel = tca_parallel,
      gene_type = gene_type,
      random_seed = random_seed,
      scale = "cpm",
      tca_version = tca_version,
      container_image = downstream_docker_image,
      docker_image = downstream_docker_image,
      cpu = manifest_cpu,
      memory = manifest_memory,
      disk_gb = manifest_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    Map[String, String] stage_images = {
      "estimation": estimation_docker_image,
      "fit": fit_docker_image,
      "export": export_docker_image,
      "downstream": downstream_docker_image
    }
    File filtered_expression = export_expression
    File? gene_type_filter_report = FilterExpressionGenes.report
    File? gene_type_filter_log = FilterExpressionGenes.log

    File? proportion_mode_validation_log = ValidateProportionMode.log

    File? estimated_proportions = MergeHspeBatches.proportions
    File? hspe_markers = PrepareHspeBatches.markers
    File? hspe_metadata = MergeHspeBatches.metadata
    File? hspe_overlap_report = PrepareHspeBatches.overlap_report
    File? transformed_lm22 = PrepareHspeBatches.transformed_lm22
    File? hspe_log = MergeHspeBatches.log
    File? hspe_sample_diagnostics = MergeHspeBatches.diagnostics

    File? proportions_lm22 = ProcessProportions.original
    File? proportions_combined = ProcessProportions.combined
    File tca_weights = CleanTcaModel.tca_weights
    File? cell_group_filter_report = ProcessProportions.filter_report
    File? proportions_log = ProcessProportions.log

    File tca_model = CleanTcaModel.model
    File? tca_model_unfiltered = FitTca.model
    File tca_numerical_excluded_genes = CleanTcaModel.excluded_genes
    File tca_cleanup_log = CleanTcaModel.log
    File? tca_model_log = FitTca.model_log
    File? tca_excluded_genes = FitTca.excluded_genes
    File? fit_tca_log = FitTca.log

    Array[File] cell_type_beds = ExportTcaBeds.cell_type_beds
    File cell_type_bed_inventory = ExportTcaBeds.cell_type_bed_inventory
    File cell_type_gene_summary = SummarizeCellTypeBeds.summary
    File gene_summary_log = SummarizeCellTypeBeds.log
    Array[File] filtered_cell_type_beds = FilterCellTypeBeds.filtered_beds
    File filtered_cell_type_bed_inventory = FilterCellTypeBeds.filtered_inventory
    File negative_expression_summary = FilterCellTypeBeds.negative_summary
    File reference_gene_comparison = FilterCellTypeBeds.gene_comparison
    File reference_filter_metrics = FilterCellTypeBeds.filter_metrics
    Array[File] reference_filter_plots = FilterCellTypeBeds.plots
    File reference_filter_log = FilterCellTypeBeds.log
    File? haemopedia_reference_summary = PrepareHaemopedia.summary
    File? haemopedia_reference_samples = PrepareHaemopedia.samples
    File? haemopedia_reference_metadata = PrepareHaemopedia.metadata
    File reconstruction_by_sample = ExportTcaBeds.reconstruction_by_sample
    File qc_summary = BuildManifest.qc_summary
    File qc_plots = ExportTcaBeds.qc_plots
    File export_log = ExportTcaBeds.log
    File export_detail_log = ExportTcaBeds.export_detail_log
    File output_manifest = BuildManifest.output_manifest
    File output_inventory = BuildManifest.provenance
    File manifest_log = BuildManifest.log
    File effective_parameters_file = BuildManifest.effective_parameters_file
  }
}
