#!/usr/bin/env Rscript

# Format preliminary methylation phenotype PCs and optional additional
# covariates for residualized local-CpG correlation analysis.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

option_list <- list(
    make_option("--PhenotypePCs", type = "character", help = "Sample-by-PC phenotype-PC TSV [required]"),
    make_option("--AdditionalCovariates", type = "character", default = NULL,
                help = "Optional sample_id-by-covariate TSV"),
    make_option("--OutputFile", type = "character", help = "TensorQTL-format covariate TSV [required]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$PhenotypePCs) || is.null(opt$OutputFile)) {
    stop("--PhenotypePCs and --OutputFile are required")
}
if (!file.exists(opt$PhenotypePCs)) stop("Phenotype-PC file does not exist: ", opt$PhenotypePCs)
if (!is.null(opt$AdditionalCovariates) && !file.exists(opt$AdditionalCovariates)) {
    stop("Additional-covariates file does not exist: ", opt$AdditionalCovariates)
}

read_sample_covariates <- function(input, label) {
    table <- if (is.character(input)) fread(input, check.names = FALSE) else copy(input)
    if (!("sample_id" %in% names(table))) {
        stop(label, " must contain a sample_id column")
    }
    table[, sample_id := as.character(sample_id)]
    if (anyNA(table$sample_id) || any(!nzchar(table$sample_id)) || anyDuplicated(table$sample_id)) {
        stop(label, " has missing, blank, or duplicate sample IDs")
    }
    if (ncol(table) < 2L) stop(label, " must contain at least one covariate column")
    value_columns <- setdiff(names(table), "sample_id")
    raw_values <- as.matrix(table[, ..value_columns])
    numeric_values <- suppressWarnings(matrix(as.numeric(raw_values), nrow = nrow(table), dimnames = list(NULL, value_columns)))
    if (any(is.na(numeric_values) & !is.na(raw_values)) || any(!is.finite(numeric_values))) {
        stop(label, " contains missing, non-numeric, or non-finite covariate values")
    }
    set(table, j = value_columns, value = as.data.table(numeric_values))
    table
}

phenotype_pcs <- fread(opt$PhenotypePCs, check.names = FALSE)
if (!("ID" %in% names(phenotype_pcs))) stop("Phenotype-PC file must contain an ID column")
setnames(phenotype_pcs, "ID", "sample_id")
phenotype_pcs <- read_sample_covariates(phenotype_pcs, "Phenotype-PC file")

if (!is.null(opt$AdditionalCovariates)) {
    additional <- read_sample_covariates(opt$AdditionalCovariates, "Additional-covariates file")
    if (!setequal(phenotype_pcs$sample_id, additional$sample_id)) {
        stop("Additional-covariates sample IDs must exactly match phenotype-PC sample IDs")
    }
    setorder(phenotype_pcs, sample_id)
    setorder(additional, sample_id)
    duplicate_columns <- intersect(setdiff(names(phenotype_pcs), "sample_id"), setdiff(names(additional), "sample_id"))
    if (length(duplicate_columns) > 0L) {
        stop("Phenotype PCs and additional covariates share column name(s): ", paste(duplicate_columns, collapse = ", "))
    }
    covariates <- cbind(phenotype_pcs, additional[, !"sample_id"])
} else {
    covariates <- phenotype_pcs
    setorder(covariates, sample_id)
}

sample_ids <- covariates$sample_id
covariate_names <- setdiff(names(covariates), "sample_id")
matrix_values <- as.matrix(covariates[, ..covariate_names])
output <- data.table(ID = covariate_names)
for (sample_index in seq_along(sample_ids)) {
    output[, (sample_ids[[sample_index]]) := matrix_values[sample_index, ]]
}

output_dir <- dirname(opt$OutputFile)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
fwrite(output, opt$OutputFile, sep = "\t", na = "NA")
message("Wrote ", length(covariate_names), " correlation covariate(s) across ", length(sample_ids), " samples: ", opt$OutputFile)
