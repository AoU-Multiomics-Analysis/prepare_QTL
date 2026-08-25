#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(RNOmni)
    library(WGCNA)
})

option_list <- list(
    make_option("--MethylationBed", type = "character",
                help = "Merged methylation BED with four interval columns followed by sample columns [required]"),
    make_option("--SampleList", type = "character",
                help = "TSV whose first column contains analysis sample IDs [required]"),
    make_option("--OutputPrefix", type = "character",
                help = "Prefix for output files [required]"),
    make_option("--MissingnessThreshold", type = "double", default = 5,
                help = "Remove features with missingness greater than or equal to this percentage [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

required_options <- c("MethylationBed", "SampleList", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--MethylationBed, --SampleList, and --OutputPrefix are required")
}
if (!is.finite(opt$MissingnessThreshold) || opt$MissingnessThreshold <= 0 || opt$MissingnessThreshold > 100) {
    stop("--MissingnessThreshold must be in (0, 100]")
}

strip_methylation_suffix <- function(ids) {
    sub("\\.meth_region_stats$", "", ids)
}

transform_feature <- function(values, method) {
    if (method == "INT") return(RNOmni::RankNorm(values))
    feature_sd <- stats::sd(values)
    if (!is.finite(feature_sd) || feature_sd == 0) return(rep(0, length(values)))
    as.numeric(scale(values, center = TRUE, scale = TRUE))
}

write_outlier_table <- function(outliers, output_file) {
    fwrite(outliers, output_file, sep = "\t", na = "NA")
}

remove_connectivity_outliers <- function(phenotype_values, output_file, transform_label) {
    outlier_file <- sub("\\.bed\\.gz$", ".connectivity_outliers.tsv", output_file)
    n_samples_before <- ncol(phenotype_values)
    empty_outliers <- data.table(SampleID = character(), Z_score = numeric())

    message("Computing sample connectivity for the ", transform_label, " matrix")
    if (n_samples_before < 3 || nrow(phenotype_values) < 2) {
        message("The ", transform_label, " matrix has too little data for connectivity analysis; keeping all ",
                n_samples_before, " samples")
        write_outlier_table(empty_outliers, outlier_file)
        return(phenotype_values)
    }

    phenotype_frame <- as.data.frame(phenotype_values, check.names = FALSE)
    adjacency <- 0.5 + 0.5 * WGCNA::bicor(phenotype_frame, use = "pairwise.complete.obs")
    adjacency[!is.finite(adjacency)] <- 0
    connectivity <- WGCNA::fundamentalNetworkConcepts(adjacency)$Connectivity
    connectivity_sd <- stats::sd(connectivity, na.rm = TRUE)

    if (!is.finite(connectivity_sd) || connectivity_sd == 0) {
        message("The ", transform_label, " connectivity scores have zero or undefined variance; keeping all ",
                n_samples_before, " samples")
        write_outlier_table(empty_outliers, outlier_file)
        return(phenotype_values)
    }

    connectivity_z <- (connectivity - mean(connectivity, na.rm = TRUE)) / connectivity_sd
    outliers <- data.table(
        SampleID = names(connectivity_z),
        Z_score = as.numeric(connectivity_z)
    )[Z_score < -3]
    setorder(outliers, Z_score, SampleID)
    write_outlier_table(outliers, outlier_file)

    kept_samples <- setdiff(colnames(phenotype_values), outliers$SampleID)
    message("Connectivity filtering for the ", transform_label, " matrix removed ", nrow(outliers),
            " of ", n_samples_before, " samples; ", length(kept_samples), " samples remain")
    phenotype_values[, kept_samples, drop = FALSE]
}

write_methylation_bed <- function(metadata, values, output_file, transform_label,
                                  transform_method = NULL, remove_outliers = TRUE) {
    if (is.null(transform_method)) {
        message("Preparing the raw methylation matrix")
        output_values <- values
    } else {
        message("Applying the ", transform_label, " transformation across samples for each feature")
        output_values <- t(vapply(
            seq_len(nrow(values)),
            function(row_index) transform_feature(as.numeric(values[row_index, ]), transform_method),
            FUN.VALUE = numeric(ncol(values))
        ))
        dimnames(output_values) <- dimnames(values)
    }

    if (remove_outliers) {
        output_values <- remove_connectivity_outliers(output_values, output_file, transform_label)
    } else {
        message("Skipping connectivity filtering for the raw matrix; keeping all ", ncol(output_values), " samples")
    }

    output <- cbind(metadata, as.data.table(output_values))
    message("Writing the ", transform_label, " BED to ", output_file)
    fwrite(output, output_file, sep = "\t", na = "NA")
}

message("Reading the sample list from ", opt$SampleList)
sample_list <- fread(opt$SampleList, header = FALSE, colClasses = "character")
if (ncol(sample_list) < 1 || nrow(sample_list) < 1) stop("--SampleList must contain at least one sample")
sample_ids <- strip_methylation_suffix(sample_list[[1]])
recognized_headers <- c("sample_id", "SampleID", "ID")
if (length(sample_ids) > 0 && sample_ids[[1]] %chin% recognized_headers) sample_ids <- sample_ids[-1]
if (length(sample_ids) < 2) stop("--SampleList must contain at least two sample IDs")
if (anyNA(sample_ids) || any(!nzchar(sample_ids))) stop("--SampleList contains an empty sample ID")
if (anyDuplicated(sample_ids)) stop("--SampleList contains duplicate sample IDs after suffix removal")
message("The sample list contains ", length(sample_ids), " samples")

message("Reading the merged methylation BED from ", opt$MethylationBed)
methylation <- fread(opt$MethylationBed, header = TRUE, na.strings = c("NA", "NaN", "."))
if (ncol(methylation) < 5) stop("--MethylationBed must contain four interval columns and at least one sample column")
setnames(methylation, names(methylation)[1:4], c("#chr", "start", "end", "phenotype_id"))
setnames(methylation, names(methylation), strip_methylation_suffix(names(methylation)))
if (anyDuplicated(names(methylation))) stop("The methylation BED contains duplicate column names after suffix removal")
if (anyNA(methylation$`#chr`) || any(!nzchar(methylation$`#chr`)) ||
    anyNA(methylation$start) || anyNA(methylation$end) ||
    anyNA(methylation$phenotype_id) || any(!nzchar(methylation$phenotype_id))) {
    stop("The methylation BED contains missing interval or phenotype metadata")
}
if (anyDuplicated(methylation$phenotype_id)) stop("The methylation BED contains duplicate phenotype IDs")

available_samples <- names(methylation)[-(1:4)]
selected_samples <- sample_ids[sample_ids %chin% available_samples]
missing_samples <- setdiff(sample_ids, available_samples)
if (length(missing_samples) > 0) {
    stop("The methylation BED is missing requested sample(s): ", paste(missing_samples, collapse = ", "))
}
if (length(selected_samples) < 2) stop("At least two sample-list IDs must match the methylation BED columns")
methylation <- methylation[, c("#chr", "start", "end", "phenotype_id", selected_samples), with = FALSE]
message("Selected ", length(selected_samples), " methylation sample columns")

for (sample_id in selected_samples) {
    converted <- suppressWarnings(as.numeric(methylation[[sample_id]]))
    introduced_missing <- is.na(converted) & !is.na(methylation[[sample_id]])
    if (any(introduced_missing)) stop("Sample column '", sample_id, "' contains a non-numeric methylation value")
    set(methylation, j = sample_id, value = converted)
}

sample_values <- as.matrix(methylation[, ..selected_samples])
percent_missing <- rowMeans(is.na(sample_values)) * 100
allowed_chromosome <- !grepl("random|X|Y|U|M", methylation$`#chr`)
keep_feature <- percent_missing < opt$MissingnessThreshold & allowed_chromosome
message("Feature filtering retained ", sum(keep_feature), " of ", nrow(methylation),
        " features using missingness < ", opt$MissingnessThreshold,
        "% and excluding random, X, Y, U, and mitochondrial contigs")
if (sum(keep_feature) < 2) {
    stop("At least two methylation features must remain after chromosome and missingness filtering for phenotype PC calculation")
}

methylation <- methylation[keep_feature]
sample_values <- sample_values[keep_feature, , drop = FALSE]
n_missing <- rowSums(is.na(sample_values))
if (nrow(sample_values) > 0 && any(n_missing > 0)) {
    for (row_index in which(n_missing > 0)) {
        feature_mean <- mean(sample_values[row_index, ], na.rm = TRUE)
        if (!is.finite(feature_mean)) stop("Cannot impute feature '", methylation$phenotype_id[[row_index]], "' because it has no observed values")
        sample_values[row_index, is.na(sample_values[row_index, ])] <- feature_mean
    }
}
message("Mean-imputed ", sum(n_missing), " value(s) across ", sum(n_missing > 0), " retained features")

metadata <- methylation[, .(`#chr`, start, end, phenotype_id)]
raw_output <- paste0(opt$OutputPrefix, ".methylation.raw.bed.gz")
int_output <- paste0(opt$OutputPrefix, ".methylation.INT.bed.gz")
scaled_output <- paste0(opt$OutputPrefix, ".methylation.scaled.bed.gz")

write_methylation_bed(metadata, sample_values, raw_output, "raw", remove_outliers = FALSE)
write_methylation_bed(metadata, sample_values, int_output, "INT", transform_method = "INT")
write_methylation_bed(metadata, sample_values, scaled_output, "scaled", transform_method = "scaled")
message("Methylation preparation completed")
