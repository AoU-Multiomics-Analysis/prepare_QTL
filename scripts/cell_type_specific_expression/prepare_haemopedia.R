#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 2L) stop("Usage: prepare_haemopedia.R COUNTS OUTPUT_DIR", call. = FALSE)
  input_path <- arguments[[1L]]
  output_dir <- arguments[[2L]]
  message(sprintf("stage=prepare_haemopedia start_time=%s", tensor_utc_time()))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  counts <- read_reference_counts(input_path)
  prepared <- prepare_reference_matrix(counts)
  summary <- summarize_reference_expression(prepared$cpm, prepared$sample_map)
  readr::write_tsv(summary, file.path(output_dir, "reference_summary.tsv.gz"))
  readr::write_tsv(prepared$samples, file.path(output_dir, "reference_samples.tsv"))
  metadata <- list(
    input = list(path = input_path, sha256 = digest::digest(file = input_path, algo = "sha256", serialize = FALSE),
                 n_genes = nrow(counts), n_samples = ncol(counts)),
    normalization = "edgeR TMM on the full raw count matrix, followed by linear CPM",
    edgeR_version = as.character(utils::packageVersion("edgeR")),
    gene_matching = "Anchored numeric Ensembl version suffixes are removed for matching only.",
    mapping_caveats = as.list(reference_mapping_caveats),
    created_utc = tensor_utc_time()
  )
  jsonlite::write_json(metadata, file.path(output_dir, "reference_metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")
  message(sprintf("stage=prepare_haemopedia genes=%d samples=%d completion_time=%s",
                  nrow(counts), ncol(counts), tensor_utc_time()))
}, error = function(error) {
  message(sprintf("stage=prepare_haemopedia status=failed utc_time=%s message=%s",
                  tensor_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})
