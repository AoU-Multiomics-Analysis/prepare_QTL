#!/usr/bin/env Rscript

# Validate a compact cohort manifest and split it into one GCS-path list per
# chromosome for workflow-internal localization.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

option_list <- list(
    make_option("--CohortManifest", type = "character", help = "TSV with sample_id, sample_qc, and autosome01..autosome22 [required]"),
    make_option("--OutputDir", type = "character", default = "manifest_lists", help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$CohortManifest) || !file.exists(opt$CohortManifest)) stop("--CohortManifest must be an existing file")

manifest <- fread(opt$CohortManifest, colClasses = "character", check.names = FALSE)
autosome_columns <- sprintf("autosome%02d", seq_len(22L))
required_columns <- c("sample_id", "sample_qc", autosome_columns)
if (!identical(names(manifest), required_columns)) {
    stop("Cohort manifest columns must be exactly: ", paste(required_columns, collapse = ", "))
}
if (nrow(manifest) < 1L) stop("Cohort manifest has no samples")
for (column in required_columns) {
    if (anyNA(manifest[[column]]) || any(!nzchar(manifest[[column]]))) stop("Cohort manifest has an empty value in ", column)
}
if (anyDuplicated(manifest$sample_id)) stop("Cohort manifest has duplicate sample_id values")

dir.create(opt$OutputDir, recursive = TRUE, showWarnings = FALSE)
fwrite(manifest[, .(sample_qc)], file.path(opt$OutputDir, "sample_qc_paths.list"), col.names = FALSE)
for (column in autosome_columns) {
    fwrite(manifest[, ..column], file.path(opt$OutputDir, paste0(column, "_paths.list")), col.names = FALSE)
}
message("Validated ", nrow(manifest), " samples and wrote chromosome path lists to ", opt$OutputDir)
