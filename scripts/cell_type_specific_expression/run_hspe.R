#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(script_path), "bootstrap.R"))

utc_time <- function() format(Sys.time(), tz = "UTC", usetz = TRUE)

run_hspe_stage <- function() {
  option_list <- list(
    optparse::make_option(
      "--expression",
      type = "character",
      help = "Coordinate-preserving BED of non-negative linear CPM values."
    ),
    optparse::make_option(
      "--log2-pseudocount",
      dest = "log2_pseudocount",
      type = "double",
      default = 0,
      help = paste(
        "Non-negative pseudocount added to bulk CPM and LM22 before",
        "their hspe log2 transforms."
      )
    ),
    optparse::make_option(
      "--gtf",
      type = "character",
      help = "GTF used to map expression gene_id values to LM22 gene symbols."
    ),
    optparse::make_option(
      "--lm22",
      type = "character",
      help = paste(
        "Standard positive linear LM22 matrix TSV with the official",
        "'Gene symbol' header or the canonical 'gene_symbol' header."
      )
    ),
    optparse::make_option(
      "--min-overlap",
      dest = "min_overlap",
      type = "double",
      default = pipeline_defaults()$min_lm22_overlap,
      help = "Minimum fraction of LM22 genes shared with bulk expression."
    ),
    optparse::make_option(
      "--marker-fraction",
      dest = "marker_fraction",
      type = "double",
      default = pipeline_defaults()$marker_fraction,
      help = "Fraction of reference genes used as hspe markers."
    ),
    optparse::make_option(
      "--marker-method",
      dest = "marker_method",
      type = "character",
      default = "ratio",
      help = "hspe marker method; only ratio is supported."
    ),
    optparse::make_option(
      "--quantile-normalize",
      dest = "quantile_normalize",
      action = "store_true",
      default = FALSE,
      help = "Jointly quantile normalize LM22 and bulk profiles for hspe only."
    ),
    optparse::make_option(
      "--random-seed",
      dest = "random_seed",
      type = "integer",
      default = 20260901L,
      help = "Random seed for HSPE optimization with DEoptimR."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for hspe outputs."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("expression", "gtf", "lm22", "output_dir")
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

  message(sprintf("stage=hspe utc_start=%s", utc_time()))
  lm22_linear <- standardize_lm22(read_lm22_matrix(options$lm22))
  expression <- read_expression_bed(options$expression)
  annotation <- read_gtf_gene_annotation(options$gtf)
  hspe_expression <- make_hspe_expression(
    expression,
    annotation,
    options$log2_pseudocount
  )
  bulk_log <- hspe_expression$log_expression
  message(sprintf(
    "stage=hspe bulk_dimensions=genes:%d samples:%d lm22_dimensions=genes:%d cell_types:%d",
    nrow(bulk_log), ncol(bulk_log), nrow(lm22_linear), ncol(lm22_linear)
  ))
  message(sprintf(
    paste0(
      "stage=hspe settings=min_overlap:%.3f marker_fraction:%.3f ",
      "marker_method:%s quantile_normalize:%s log2_pseudocount:%g random_seed:%d optimizer:DEoptimR"
    ),
    options$min_overlap,
    options$marker_fraction,
    options$marker_method,
    options$quantile_normalize,
    hspe_expression$log2_pseudocount,
    options$random_seed
  ))
  inputs <- prepare_hspe_inputs(
    bulk_log = bulk_log,
    lm22_linear = lm22_linear,
    min_overlap = options$min_overlap,
    quantile_normalize = options$quantile_normalize,
    log2_pseudocount = hspe_expression$log2_pseudocount
  )
  message(sprintf(
    "stage=hspe overlap=shared_genes:%d overlap_fraction:%.3f",
    inputs$overlap_count,
    inputs$overlap_fraction
  ))
  fit <- estimate_hspe(
    inputs,
    marker_fraction = options$marker_fraction,
    marker_method = options$marker_method,
    random_seed = options$random_seed
  )
  fit$metadata$log2_pseudocount <- hspe_expression$log2_pseudocount
  marker_counts <- unlist(fit$metadata$marker_counts, use.names = TRUE)
  message(sprintf("stage=hspe package=hspe version=%s", fit$metadata$hspe_version))
  message(sprintf(
    "stage=hspe markers_per_cell_type=%s",
    paste(sprintf("%s:%d", names(marker_counts), marker_counts), collapse = ",")
  ))

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    proportions = file.path(options$output_dir, "hspe_proportions.tsv"),
    markers = file.path(options$output_dir, "hspe_markers.tsv"),
    metadata = file.path(options$output_dir, "hspe_metadata.json"),
    overlap = file.path(options$output_dir, "hspe_overlap.tsv"),
    lm22_log = file.path(options$output_dir, "hspe_lm22_log.tsv.gz")
  )
  message(sprintf(
    "stage=hspe output_paths=%s",
    paste(unlist(output_paths, use.names = FALSE), collapse = ",")
  ))
  write_numeric_matrix(fit$proportions, output_paths$proportions, "sample_id")
  readr::write_tsv(fit$markers, output_paths$markers, na = "")
  jsonlite::write_json(fit$metadata, output_paths$metadata, auto_unbox = TRUE, pretty = TRUE)
  readr::write_tsv(inputs$overlap_report, output_paths$overlap, na = "")
  write_numeric_matrix(inputs$transformed_lm22, output_paths$lm22_log, "gene_symbol")
  message(sprintf("stage=hspe utc_complete=%s", utc_time()))
}

tryCatch(
  run_hspe_stage(),
  error = function(error) {
    message(sprintf(
      "stage=hspe status=failed utc_time=%s message=%s",
      utc_time(), conditionMessage(error)
    ))
    quit(status = 1L)
  }
)
