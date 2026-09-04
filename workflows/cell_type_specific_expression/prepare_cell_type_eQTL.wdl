version 1.0

import "deconvolution.wdl" as deconvolution
import "tasks/integration.wdl" as integration
import "../expression/prepare_eQTL.wdl" as eqtl

workflow PrepareCellTypeEqtlWorkflow {
  input {
    File expression
    File gtf
    File lm22
    File? precomputed_proportions
    File? deconvolution_covariates
    File SampleList
    File AdditionalCovariates
    String OutputPrefix

    String deconvolution_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main"
    String qtl_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
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
    Int scatter_cpu = 1
    String scatter_memory = "4 GB"
    Int scatter_disk_gb = 20
    Int eqtl_cpu = 8
    Int eqtl_memory = 64
    Int eqtl_disk_gb = 200
  }

  call deconvolution.CellTypeDeconvolution as CellTypeDeconvolution {
    input:
      expression = expression,
      gtf = gtf,
      lm22 = lm22,
      precomputed_proportions = precomputed_proportions,
      covariates = deconvolution_covariates,
      deconvolution_docker_image = deconvolution_docker_image,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries,
      min_lm22_overlap = min_lm22_overlap,
      hspe_marker_fraction = hspe_marker_fraction,
      hspe_quantile_normalize = hspe_quantile_normalize,
      hspe_batch_size = hspe_batch_size,
      hspe_batch_memory = hspe_batch_memory,
      hspe_batch_disk_gb = hspe_batch_disk_gb,
      group_mean_threshold = group_mean_threshold,
      zero_floor = zero_floor,
      tca_max_iters = tca_max_iters,
      tca_parallel = tca_parallel,
      random_seed = random_seed,
      log2_pseudocount = log2_pseudocount,
      gene_type = gene_type,
      hspe_cpu = hspe_cpu,
      hspe_memory = hspe_memory,
      hspe_disk_gb = hspe_disk_gb,
      proportions_cpu = proportions_cpu,
      proportions_memory = proportions_memory,
      proportions_disk_gb = proportions_disk_gb,
      fit_cpu = fit_cpu,
      fit_memory = fit_memory,
      fit_disk_gb = fit_disk_gb,
      export_cpu = export_cpu,
      export_memory = export_memory,
      export_disk_gb = export_disk_gb,
      gene_summary_memory = gene_summary_memory,
      manifest_cpu = manifest_cpu,
      manifest_memory = manifest_memory,
      manifest_disk_gb = manifest_disk_gb
  }

  call integration.PrepareScatterInputs as PrepareScatterInputs {
    input:
      cell_type_bed_inventory = CellTypeDeconvolution.cell_type_bed_inventory,
      cell_type_beds = CellTypeDeconvolution.cell_type_beds,
      output_prefix = OutputPrefix,
      docker_image = deconvolution_docker_image,
      cpu = scatter_cpu,
      memory = scatter_memory,
      disk_gb = scatter_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  scatter (index in range(length(PrepareScatterInputs.cell_types))) {
    call eqtl.eQTLPrepareData as PrepareCellTypeEqtl {
      input:
        OutputPrefix = PrepareScatterInputs.output_prefixes[index],
        CpmBed = PrepareScatterInputs.expression_beds[index],
        SampleList = SampleList,
        AdditionalCovariates = AdditionalCovariates,
        ResidualizeNormalizedInputs = false,
        DockerImage = qtl_docker_image,
        preemptible_attempts = preemptible_attempts,
        max_retries = max_retries,
        memory = eqtl_memory,
        disk_space = eqtl_disk_gb,
        num_threads = eqtl_cpu
    }

    File int_merged_covariate = select_first([PrepareCellTypeEqtl.IntQtlCovariates])
    File scaled_merged_covariate = select_first([PrepareCellTypeEqtl.ScaledQtlCovariates])
  }

  call integration.BuildQtlManifest as BuildQtlManifest {
    input:
      cell_types = PrepareScatterInputs.cell_types,
      cell_type_slugs = PrepareScatterInputs.cell_type_slugs,
      int_beds = PrepareCellTypeEqtl.IntBedFile,
      scaled_beds = PrepareCellTypeEqtl.ScaledBedFile,
      int_phenotype_pcs = PrepareCellTypeEqtl.IntPhenotypePCsOut,
      int_phenotype_pcs_all = PrepareCellTypeEqtl.IntPhenotypePCsAllOut,
      scaled_phenotype_pcs = PrepareCellTypeEqtl.ScaledPhenotypePCsOut,
      scaled_phenotype_pcs_all = PrepareCellTypeEqtl.ScaledPhenotypePCsAllOut,
      int_merged_covariates = int_merged_covariate,
      scaled_merged_covariates = scaled_merged_covariate,
      int_connectivity_outliers = PrepareCellTypeEqtl.IntConnectivityOutliers,
      scaled_connectivity_outliers = PrepareCellTypeEqtl.ScaledConnectivityOutliers,
      docker_image = deconvolution_docker_image,
      cpu = scatter_cpu,
      memory = scatter_memory,
      disk_gb = scatter_disk_gb,
      preemptible_attempts = preemptible_attempts,
      max_retries = max_retries
  }

  output {
    File filtered_expression = CellTypeDeconvolution.filtered_expression
    File gene_type_filter_report = CellTypeDeconvolution.gene_type_filter_report
    File gene_type_filter_log = CellTypeDeconvolution.gene_type_filter_log

    File proportion_mode_validation_log = CellTypeDeconvolution.proportion_mode_validation_log

    File? estimated_proportions = CellTypeDeconvolution.estimated_proportions
    File? hspe_markers = CellTypeDeconvolution.hspe_markers
    File? hspe_metadata = CellTypeDeconvolution.hspe_metadata
    File? hspe_overlap_report = CellTypeDeconvolution.hspe_overlap_report
    File? transformed_lm22 = CellTypeDeconvolution.transformed_lm22
    File? hspe_log = CellTypeDeconvolution.hspe_log
    File? hspe_sample_diagnostics = CellTypeDeconvolution.hspe_sample_diagnostics

    File proportions_lm22 = CellTypeDeconvolution.proportions_lm22
    File proportions_combined = CellTypeDeconvolution.proportions_combined
    File tca_weights = CellTypeDeconvolution.tca_weights
    File cell_group_filter_report = CellTypeDeconvolution.cell_group_filter_report
    File proportions_log = CellTypeDeconvolution.proportions_log

    File tca_model = CellTypeDeconvolution.tca_model
    File tca_model_unfiltered = CellTypeDeconvolution.tca_model_unfiltered
    File tca_numerical_excluded_genes = CellTypeDeconvolution.tca_numerical_excluded_genes
    File tca_cleanup_log = CellTypeDeconvolution.tca_cleanup_log
    File tca_model_log = CellTypeDeconvolution.tca_model_log
    File tca_excluded_genes = CellTypeDeconvolution.tca_excluded_genes
    File fit_tca_log = CellTypeDeconvolution.fit_tca_log

    Array[File] cell_type_beds = CellTypeDeconvolution.cell_type_beds
    File cell_type_bed_inventory = CellTypeDeconvolution.cell_type_bed_inventory
    File cell_type_gene_summary = CellTypeDeconvolution.cell_type_gene_summary
    File gene_summary_log = CellTypeDeconvolution.gene_summary_log
    File reconstruction_by_sample = CellTypeDeconvolution.reconstruction_by_sample
    File qc_summary = CellTypeDeconvolution.qc_summary
    File qc_plots = CellTypeDeconvolution.qc_plots
    File export_log = CellTypeDeconvolution.export_log
    File export_detail_log = CellTypeDeconvolution.export_detail_log
    File output_manifest = CellTypeDeconvolution.output_manifest
    File output_inventory = CellTypeDeconvolution.output_inventory
    File manifest_log = CellTypeDeconvolution.manifest_log
    File effective_parameters_file = CellTypeDeconvolution.effective_parameters_file

    File cell_type_qtl_manifest = BuildQtlManifest.manifest
    File cell_type_qtl_manifest_log = BuildQtlManifest.log
    Array[String] cell_types = PrepareScatterInputs.cell_types
    Array[String] cell_type_slugs = PrepareScatterInputs.cell_type_slugs
    Array[File] int_beds = PrepareCellTypeEqtl.IntBedFile
    Array[File] scaled_beds = PrepareCellTypeEqtl.ScaledBedFile
    Array[File] int_phenotype_pcs = PrepareCellTypeEqtl.IntPhenotypePCsOut
    Array[File] int_phenotype_pcs_all = PrepareCellTypeEqtl.IntPhenotypePCsAllOut
    Array[File] scaled_phenotype_pcs = PrepareCellTypeEqtl.ScaledPhenotypePCsOut
    Array[File] scaled_phenotype_pcs_all = PrepareCellTypeEqtl.ScaledPhenotypePCsAllOut
    Array[File] int_merged_covariates = int_merged_covariate
    Array[File] scaled_merged_covariates = scaled_merged_covariate
    Array[File] int_connectivity_outliers = PrepareCellTypeEqtl.IntConnectivityOutliers
    Array[File] scaled_connectivity_outliers = PrepareCellTypeEqtl.ScaledConnectivityOutliers
  }
}
