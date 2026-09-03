#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

run_filter_expression_genes <- function() {
  option_list <- list(
    optparse::make_option(
      "--expression",
      type = "character",
      help = "Coordinate-preserving BED of non-negative linear CPM values."
    ),
    optparse::make_option(
      "--gtf",
      type = "character",
      help = "GTF with gene_id and gene_type or gene_biotype attributes."
    ),
    optparse::make_option(
      "--gene-types",
      dest = "gene_types",
      type = "character",
      help = "Comma-separated GTF gene types to retain."
    ),
    optparse::make_option(
      "--log2-pseudocount",
      dest = "log2_pseudocount",
      type = "double",
      default = 0,
      help = "Non-negative pseudocount used by downstream log2 transforms."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for the filtered BED and gene filter report."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("expression", "gtf", "gene_types", "output_dir")
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf("Missing required options: %s", paste(missing_options, collapse = ", ")),
      call. = FALSE
    )
  }

  gene_types <- parse_gene_types(options$gene_types)
  message(sprintf(
    "stage=filter_expression_genes utc_start=%s gene_types=%s",
    utc_time(),
    paste(gene_types, collapse = ",")
  ))
  expression <- read_expression_bed(
    options$expression,
    options$log2_pseudocount
  )
  annotation <- read_gtf_gene_annotation(
    options$gtf,
    require_gene_name = FALSE
  )
  filtered <- filter_expression_by_gene_types(
    expression,
    annotation,
    gene_types
  )

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  expression_path <- file.path(options$output_dir, "filtered_expression.bed.gz")
  report_path <- file.path(options$output_dir, "gene_type_filter_report.tsv")
  write_expression_bed(
    expression_path,
    filtered$expression$coordinates,
    filtered$expression$cpm
  )
  readr::write_tsv(filtered$report, report_path, na = "")
  message(sprintf(
    paste0(
      "stage=filter_expression_genes input_genes=%d retained_genes=%d ",
      "output_paths=%s,%s completion_time=%s"
    ),
    nrow(expression$coordinates),
    nrow(filtered$expression$coordinates),
    expression_path,
    report_path,
    utc_time()
  ))
}

tryCatch(
  run_filter_expression_genes(),
  error = function(error) {
    message(sprintf(
      "stage=filter_expression_genes status=failed utc_time=%s message=%s",
      utc_time(),
      conditionMessage(error)
    ))
    quit(save = "no", status = 1L)
  }
)
