#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--expression", type = "character"),
    optparse::make_option("--gtf", type = "character"),
    optparse::make_option("--lm22", type = "character"),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character"),
    optparse::make_option("--batch-size", dest = "batch_size", type = "integer", default = 100L),
    optparse::make_option("--random-seed", dest = "random_seed", type = "integer", default = 20260901L),
    optparse::make_option("--marker-fraction", dest = "marker_fraction", type = "double", default = .1),
    optparse::make_option("--min-overlap", dest = "min_overlap", type = "double", default = .8),
    optparse::make_option("--log2-pseudocount", dest = "log2_pseudocount", type = "double", default = 0),
    optparse::make_option("--quantile-normalize", dest = "quantile_normalize",
                          action = "store_true", default = FALSE)
  )))
  required <- c("expression", "gtf", "lm22", "output_dir")
  if (any(purrr::map_lgl(options[required], ~ is.null(.x) || !nzchar(.x)))) {
    stop("expression, gtf, lm22, and output-dir are required", call. = FALSE)
  }
  message(sprintf("stage=prepare_hspe_batches input_paths=%s,%s,%s",
                  options$expression, options$gtf, options$lm22))
  expression <- read_expression_bed(options$expression)
  annotation <- read_gtf_gene_annotation(options$gtf)
  bulk <- make_hspe_expression(expression, annotation, options$log2_pseudocount)
  inputs <- prepare_hspe_inputs(bulk$log_expression, read_lm22_matrix(options$lm22),
    min_overlap = options$min_overlap, quantile_normalize = options$quantile_normalize,
    log2_pseudocount = options$log2_pseudocount)
  prepared <- prepare_hspe_batches(inputs, options$batch_size,
                                   options$marker_fraction, options$random_seed)
  prepared$shared$metadata$log2_pseudocount <- options$log2_pseudocount
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(prepared$shared, file.path(options$output_dir, "hspe_prepared.rds"))
  purrr::iwalk(prepared$batches, function(batch, index) {
    saveRDS(batch, file.path(options$output_dir, sprintf("batch_%06d.rds", index)))
  })
  readr::write_tsv(prepared$shared$marker_table, file.path(options$output_dir, "hspe_markers.tsv"))
  readr::write_tsv(inputs$overlap_report, file.path(options$output_dir, "hspe_overlap.tsv"))
  write_numeric_matrix(inputs$transformed_lm22,
                        file.path(options$output_dir, "hspe_lm22_log.tsv.gz"), "gene_symbol")
  message(sprintf("stage=prepare_hspe_batches dimensions=samples:%d,markers:%d,batches:%d outputs=%s",
    nrow(inputs$Y), ncol(prepared$shared$references), length(prepared$batches), options$output_dir))
}, error = function(error) {
  message(sprintf("stage=prepare_hspe_batches status=failed message=%s", conditionMessage(error)))
  quit(status = 1L)
})
