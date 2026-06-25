library(tidyverse)
library(data.table)
library(optparse)


######## FUNCTIONS ########

is_present <- function(x){
    !is.null(x) && !is.na(x) && nzchar(x) && !tolower(x) %in% c("none", "null", "na")
}

default_output_file <- function(input_bed){
    if (str_detect(input_bed, "\\.bed\\.gz$")) {
        return(str_replace(input_bed, "\\.bed\\.gz$", ".residualized.bed.gz"))
    }
    paste0(input_bed, ".residualized.bed.gz")
}

check_required_bed <- function(bed_df){
    if (ncol(bed_df) < 5) {
        stop("InputBed must contain at least 4 metadata columns and one sample column")
    }
}

to_numeric_matrix <- function(df, matrix_name){
    raw_values <- as.matrix(df)
    numeric_values <- suppressWarnings(
        matrix(as.numeric(raw_values), nrow = nrow(raw_values), dimnames = dimnames(raw_values))
    )
    non_numeric <- is.na(numeric_values) & !is.na(raw_values)
    if (any(non_numeric)) {
        stop(paste0(matrix_name, " contains non-numeric values"))
    }
    numeric_values
}

scale_rows <- function(mat){
    row_centers <- rowMeans(mat, na.rm = TRUE)
    row_centers[is.nan(row_centers)] <- NA_real_

    centered <- sweep(mat, 1, row_centers, FUN = "-")
    row_sds <- apply(centered, 1, sd, na.rm = TRUE)
    zero_variance <- is.na(row_sds) | row_sds == 0

    if (any(zero_variance)) {
        message(paste0(sum(zero_variance), " phenotype rows have zero or undefined variance after residualization/scaling"))
    }

    row_sds[zero_variance] <- NA_real_
    sweep(centered, 1, row_sds, FUN = "/")
}

load_covariates <- function(covariate_file){
    covariates <- fread(
        covariate_file,
        data.table = FALSE,
        check.names = FALSE,
        na.strings = c("", "NA", "NaN", ".")
    )

    if (ncol(covariates) < 2) {
        stop("Covariates file must contain one covariate ID column and at least one sample column")
    }

    covariate_names <- as.character(covariates[[1]])
    if (any(is.na(covariate_names)) || any(!nzchar(covariate_names))) {
        stop("Covariate ID column contains missing or blank covariate names")
    }

    covariate_values <- covariates[, -1, drop = FALSE]
    rownames(covariate_values) <- make.unique(covariate_names)

    covariate_df <- as.data.frame(t(as.matrix(covariate_values)), stringsAsFactors = FALSE)
    colnames(covariate_df) <- make.names(rownames(covariate_values), unique = TRUE)
    covariate_df <- type.convert(covariate_df, as.is = TRUE)

    covariate_df
}

build_design_matrix <- function(covariate_df){
    complete_samples <- complete.cases(covariate_df)
    if (!all(complete_samples)) {
        message(paste0(sum(!complete_samples), " samples removed because covariates contain missing values"))
    }

    covariate_df <- covariate_df[complete_samples, , drop = FALSE]
    design_matrix <- model.matrix(~ ., data = covariate_df)

    design_rank <- qr(design_matrix)$rank
    if (nrow(design_matrix) <= design_rank) {
        stop("Not enough samples to residualize: sample count must exceed covariate design rank")
    }

    list(
        covariates = covariate_df,
        design = design_matrix
    )
}

residualize_matrix <- function(phenotype_matrix, design_matrix){
    if (!anyNA(phenotype_matrix)) {
        fit <- lm.fit(x = design_matrix, y = t(phenotype_matrix))
        residuals <- t(fit$residuals)
        rownames(residuals) <- rownames(phenotype_matrix)
        colnames(residuals) <- colnames(phenotype_matrix)
        return(residuals)
    }

    message("Phenotype matrix contains missing values; residualizing rows individually")
    residuals <- matrix(
        NA_real_,
        nrow = nrow(phenotype_matrix),
        ncol = ncol(phenotype_matrix),
        dimnames = dimnames(phenotype_matrix)
    )

    skipped_rows <- 0
    for (row_index in seq_len(nrow(phenotype_matrix))) {
        y <- phenotype_matrix[row_index, ]
        complete_y <- !is.na(y)
        row_design <- design_matrix[complete_y, , drop = FALSE]
        row_rank <- qr(row_design)$rank

        if (sum(complete_y) <= row_rank) {
            skipped_rows <- skipped_rows + 1
            next
        }

        fit <- lm.fit(x = row_design, y = y[complete_y])
        residuals[row_index, complete_y] <- fit$residuals
    }

    if (skipped_rows > 0) {
        message(paste0(skipped_rows, " phenotype rows skipped because they did not have enough complete samples"))
    }

    residuals
}


######## COMMAND LINE ARGUMENTS ########

option_list <- list(
    optparse::make_option(c("--InputBed"), type = "character", default = NULL,
                        help = "Normalized molecular phenotype BED file", metavar = "type"),
    optparse::make_option(c("--Covariates"), type = "character", default = NULL,
                        help = "Optional merged covariates TSV in tensorQTL format", metavar = "type"),
    optparse::make_option(c("--OutputFile"), type = "character", default = NULL,
                        help = "Output BED file. Defaults to inserting .residualized before .bed.gz", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

if (!is_present(opt$InputBed)) {
    stop("--InputBed is required")
}

OutputFile <- opt$OutputFile
if (!is_present(OutputFile)) {
    OutputFile <- default_output_file(opt$InputBed)
}


######## LOAD BED ########

message(paste0("Reading input BED: ", opt$InputBed))
bed_df <- fread(
    opt$InputBed,
    data.table = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN", ".")
)
check_required_bed(bed_df)

metadata_df <- bed_df[, 1:4, drop = FALSE]
sample_ids <- colnames(bed_df)[5:ncol(bed_df)]
phenotype_matrix <- to_numeric_matrix(bed_df[, sample_ids, drop = FALSE], "InputBed phenotype matrix")
rownames(phenotype_matrix) <- seq_len(nrow(phenotype_matrix))


######## RESIDUALIZE AND SCALE ########

if (is_present(opt$Covariates)) {
    message(paste0("Reading merged covariates: ", opt$Covariates))
    covariate_df <- load_covariates(opt$Covariates)

    overlapping_samples <- sample_ids[sample_ids %in% rownames(covariate_df)]
    if (length(overlapping_samples) == 0) {
        stop("No overlapping sample IDs found between InputBed and Covariates")
    }

    dropped_bed_samples <- setdiff(sample_ids, overlapping_samples)
    if (length(dropped_bed_samples) > 0) {
        message(paste0(length(dropped_bed_samples), " BED samples dropped because they were not present in covariates"))
    }

    covariate_df <- covariate_df[overlapping_samples, , drop = FALSE]
    design <- build_design_matrix(covariate_df)

    kept_samples <- rownames(design$covariates)
    phenotype_matrix <- phenotype_matrix[, kept_samples, drop = FALSE]

    message("Residualizing phenotypes onto covariates")
    residual_matrix <- residualize_matrix(phenotype_matrix, design$design)

    message("Scaling residuals")
    output_matrix <- scale_rows(residual_matrix)
} else {
    message("No covariates supplied; scaling input phenotypes without residualization")
    kept_samples <- sample_ids
    output_matrix <- scale_rows(phenotype_matrix)
}


######## WRITE OUTPUT ########

output_df <- bind_cols(
    metadata_df,
    as.data.frame(output_matrix, check.names = FALSE)
)
colnames(output_df) <- c(colnames(metadata_df), kept_samples)

message(paste0("Writing residualized BED: ", OutputFile))
fwrite(output_df, OutputFile, sep = "\t", quote = FALSE, na = "NA")
