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
    Float dtangle_marker_fraction = 0.10
    Boolean dtangle_quantile_normalize = false
    Float group_mean_threshold = 0.0001
    Float zero_floor = 0.000001
    Int tca_max_iters = 10
    Boolean tca_parallel = false
    Int random_seed = 20260901
    Float log2_pseudocount = 0.0
    Array[String] gene_type = ["protein_coding", "lncRNA"]

    Int dtangle_cpu = 4
    String dtangle_memory = "32 GB"
    Int dtangle_disk_gb = 100
    Int proportions_cpu = 2
    String proportions_memory = "16 GB"
    Int proportions_disk_gb = 50
    Int fit_cpu = 16
    String fit_memory = "256 GB"
    Int fit_disk_gb = 750
    Int export_cpu = 8
    String export_memory = "256 GB"
    Int export_disk_gb = 500
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
      dtangle_marker_fraction = dtangle_marker_fraction,
      dtangle_quantile_normalize = dtangle_quantile_normalize,
      group_mean_threshold = group_mean_threshold,
      zero_floor = zero_floor,
      tca_max_iters = tca_max_iters,
      tca_parallel = tca_parallel,
      random_seed = random_seed,
      log2_pseudocount = log2_pseudocount,
      gene_type = gene_type,
      dtangle_cpu = dtangle_cpu,
      dtangle_memory = dtangle_memory,
      dtangle_disk_gb = dtangle_disk_gb,
      proportions_cpu = proportions_cpu,
      proportions_memory = proportions_memory,
      proportions_disk_gb = proportions_disk_gb,
      fit_cpu = fit_cpu,
      fit_memory = fit_memory,
      fit_disk_gb = fit_disk_gb,
      export_cpu = export_cpu,
      export_memory = export_memory,
      export_disk_gb = export_disk_gb,
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
        Log2CpmBed = PrepareScatterInputs.expression_beds[index],
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
    File? dtangle_markers = CellTypeDeconvolution.dtangle_markers
    File? dtangle_metadata = CellTypeDeconvolution.dtangle_metadata
    File? dtangle_overlap_report = CellTypeDeconvolution.dtangle_overlap_report
    File? transformed_lm22 = CellTypeDeconvolution.transformed_lm22
    File? dtangle_log = CellTypeDeconvolution.dtangle_log

    File proportions_lm22 = CellTypeDeconvolution.proportions_lm22
    File proportions_combined = CellTypeDeconvolution.proportions_combined
    File tca_weights = CellTypeDeconvolution.tca_weights
    File cell_group_filter_report = CellTypeDeconvolution.cell_group_filter_report
    File proportions_log = CellTypeDeconvolution.proportions_log

    File tca_model = CellTypeDeconvolution.tca_model
    File tca_model_log = CellTypeDeconvolution.tca_model_log
    File tca_expression = CellTypeDeconvolution.tca_expression
    File tca_excluded_genes = CellTypeDeconvolution.tca_excluded_genes
    File fit_tca_log = CellTypeDeconvolution.fit_tca_log

    Array[File] cell_type_beds = CellTypeDeconvolution.cell_type_beds
    File cell_type_bed_inventory = CellTypeDeconvolution.cell_type_bed_inventory
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
