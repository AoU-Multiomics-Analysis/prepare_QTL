#!/usr/bin/env Rscript
# Legacy stage filename retained for the task-image release registry.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(dirname(normalizePath(script_path))), "bootstrap.R"))
tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--reference", type = "character"),
    optparse::make_option("--gtf", type = "character"),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character")
  )))
  for (key in c("reference", "gtf", "output_dir")) {
    if (is.null(options[[key]]) || !nzchar(options[[key]])) stop(paste("Missing option", key))
  }
  for (path in c(options$reference, options$gtf)) {
    if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path) || !file.exists(path) || file.access(path, 4) != 0) {
      stop(paste("Input localization error: unreadable local file", path))
    }
  }
  message(sprintf("stage=prepare_tabula_sapiens start_time=%s", tensor_utc_time()))
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  profiles <- readr::read_tsv(options$reference, show_col_types = FALSE, name_repair = "minimal")
  annotation <- read_gtf_gene_annotation(options$gtf)
  prepared <- prepare_tabula_reference(profiles, annotation)
  readr::write_tsv(prepared$summary, file.path(options$output_dir, "reference_summary.tsv.gz"))
  readr::write_tsv(prepared$samples, file.path(options$output_dir, "reference_samples.tsv"))
  metadata <- list(
    reference = options$reference,
    reference_sha256 = digest::digest(file = options$reference, algo = "sha256", serialize = FALSE),
    gtf_sha256 = digest::digest(file = options$gtf, algo = "sha256", serialize = FALSE),
    normalization = "None: supplied ComBat-corrected linear CPM; log2(1 + CPM) for comparison",
    gene_matching = "Unambiguous GTF gene_name to version-stripped gene_id; no gene aggregation",
    excluded_symbols = prepared$excluded_symbols,
    n_reference_profiles_per_cell_type = 1L,
    evidence = "Reference corrected using cohort bulk data; agreement is not independent validation",
    created_utc = tensor_utc_time())
  jsonlite::write_json(metadata, file.path(options$output_dir, "reference_metadata.json"),
    auto_unbox = TRUE, pretty = TRUE)
  message(sprintf("stage=prepare_tabula_sapiens status=complete genes=%d", dplyr::n_distinct(prepared$summary$gene_id)))
}, error = function(error) {
  message(sprintf("stage=prepare_tabula_sapiens status=failed message=%s", conditionMessage(error)))
  quit(status = 1L)
})
