#!/usr/bin/env Rscript

# Calculate full sample connectivity from locally de-correlated CpG
# representatives, then retain the passing samples in every final phenotype file.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(WGCNA)
})

open_input_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        if (nzchar(Sys.which("bgzip"))) return(pipe(sprintf("bgzip -c -d %s", shQuote(path)), open = "r"))
        return(gzfile(path, open = "rt"))
    }
    file(path, open = "rt")
}

open_output_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) return(gzfile(path, open = "wt"))
    file(path, open = "wt")
}

read_file_list <- function(path, label) {
    if (!file.exists(path)) stop(label, " list does not exist: ", path)
    files <- scan(path, what = character(), quiet = TRUE)
    if (length(files) == 0L) stop(label, " list is empty")
    list_dir <- dirname(normalizePath(path))
    relative <- !grepl("^(/|~)", files)
    files[relative] <- file.path(list_dir, files[relative])
    files <- path.expand(files)
    missing <- files[!file.exists(files)]
    if (length(missing) > 0L) stop(label, " file(s) do not exist: ", paste(missing, collapse = ", "))
    files
}

read_bed_header <- function(path) {
    connection <- open_input_connection(path)
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header) != 1L) stop("BED is empty: ", path)
    columns <- strsplit(header, "\t", fixed = TRUE)[[1]]
    if (length(columns) < 5L) stop("BED must contain four metadata columns and sample columns: ", path)
    columns
}

copy_or_filter_wide_bed <- function(input_path, output_path, keep_samples) {
    input_connection <- open_input_connection(input_path)
    on.exit(close(input_connection), add = TRUE)
    header <- readLines(input_connection, n = 1L, warn = FALSE)
    if (length(header) != 1L) stop("BED is empty: ", input_path)
    columns <- strsplit(header, "\t", fixed = TRUE)[[1]]
    sample_columns <- columns[-seq_len(4L)]
    missing_samples <- setdiff(keep_samples, sample_columns)
    if (length(missing_samples) > 0L) stop("BED is missing retained sample(s): ", paste(missing_samples, collapse = ", "))
    keep_columns <- c(seq_len(4L), match(keep_samples, columns))
    output_connection <- open_output_connection(output_path)
    on.exit(close(output_connection), add = TRUE)
    writeLines(paste(columns[keep_columns], collapse = "\t"), output_connection)
    repeat {
        lines <- readLines(input_connection, n = 250L, warn = FALSE)
        if (length(lines) == 0L) break
        filtered <- vapply(strsplit(lines, "\t", fixed = TRUE), function(fields) {
            if (length(fields) != length(columns)) stop("BED row has a different number of columns than its header: ", input_path)
            paste(fields[keep_columns], collapse = "\t")
        }, character(1))
        writeLines(filtered, output_connection)
    }
    invisible(NULL)
}

copy_or_filter_long_calls <- function(input_path, output_path, keep_samples) {
    input_connection <- open_input_connection(input_path)
    on.exit(close(input_connection), add = TRUE)
    header <- readLines(input_connection, n = 1L, warn = FALSE)
    if (length(header) != 1L) stop("Filtered-call table is empty: ", input_path)
    if (!identical(strsplit(header, "\t", fixed = TRUE)[[1]][[1]], "sample_id")) {
        stop("Filtered-call table must start with sample_id: ", input_path)
    }
    output_connection <- open_output_connection(output_path)
    on.exit(close(output_connection), add = TRUE)
    writeLines(header, output_connection)
    repeat {
        lines <- readLines(input_connection, n = 10000L, warn = FALSE)
        if (length(lines) == 0L) break
        sample_ids <- sub("\t.*$", "", lines)
        writeLines(lines[sample_ids %in% keep_samples], output_connection)
    }
    invisible(NULL)
}

load_representative_values <- function(int_bed_paths, representatives, chunk_rows) {
    reference_header <- read_bed_header(int_bed_paths[[1]])
    sample_ids <- reference_header[-seq_len(4L)]
    if (anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
        stop("INT BED contains missing, blank, or duplicate sample IDs")
    }
    values <- matrix(NA_real_, nrow = nrow(representatives), ncol = length(sample_ids),
                     dimnames = list(representatives$phenotype_id, sample_ids))
    found <- rep(FALSE, nrow(representatives))
    representative_index <- setNames(seq_len(nrow(representatives)), representatives$phenotype_id)

    for (input_path in int_bed_paths) {
        header <- read_bed_header(input_path)
        if (!identical(header, reference_header)) stop("Chromosome INT BEDs do not have identical headers: ", input_path)
        connection <- open_input_connection(input_path)
        on.exit(close(connection), add = TRUE)
        ignored_header <- readLines(connection, n = 1L, warn = FALSE)
        message("Reading representative CpGs from ", input_path)
        repeat {
            lines <- readLines(connection, n = chunk_rows, warn = FALSE)
            if (length(lines) == 0L) break
            block <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE, sep = "\t",
                           col.names = reference_header, check.names = FALSE, showProgress = FALSE)
            matching_rows <- unname(representative_index[as.character(block[[reference_header[[4]]]])])
            keep <- which(!is.na(matching_rows))
            if (length(keep) == 0L) next
            if (any(found[matching_rows[keep]])) stop("A representative phenotype ID occurs more than once in the INT BEDs")
            raw_values <- as.matrix(block[keep, ..sample_ids])
            numeric_values <- suppressWarnings(matrix(as.numeric(raw_values), nrow = length(keep), ncol = length(sample_ids)))
            if (any(is.na(numeric_values) & !is.na(raw_values)) || any(!is.finite(numeric_values))) {
                stop("Representative INT CpGs contain missing, non-numeric, or non-finite values")
            }
            values[matching_rows[keep], ] <- numeric_values
            found[matching_rows[keep]] <- TRUE
        }
        close(connection)
        on.exit(NULL, add = FALSE)
    }
    if (any(!found)) stop("Representative CpG(s) were not found in the chromosome INT BEDs: ",
                           paste(representatives$phenotype_id[!found], collapse = ", "))
    list(values = values, sample_ids = sample_ids)
}

option_list <- list(
    make_option("--IntBedList", type = "character", help = "One per-chromosome pre-connectivity INT BED path per line [required]"),
    make_option("--RepresentativeList", type = "character", help = "One per-chromosome representative-CpG TSV path per line [required]"),
    make_option("--FilteredCalls", type = "character", help = "Pre-connectivity long call table [required]"),
    make_option("--RawBed", type = "character", help = "Pre-connectivity raw methylation BED [required]"),
    make_option("--IntBed", type = "character", help = "Pre-connectivity INT methylation BED [required]"),
    make_option("--SampleQC", type = "character", help = "Pre-connectivity sample-QC TSV [required]"),
    make_option("--OutputPrefix", type = "character", help = "Output prefix [required]"),
    make_option("--ConnectivityZThreshold", type = "double", default = -3,
                help = "Samples below this connectivity Z score are removed [default: %default]"),
    make_option("--ChunkRows", type = "integer", default = 100,
                help = "INT BED rows parsed at once while extracting representatives [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("IntBedList", "RepresentativeList", "FilteredCalls", "RawBed", "IntBed", "SampleQC", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--IntBedList, --RepresentativeList, --FilteredCalls, --RawBed, --IntBed, --SampleQC, and --OutputPrefix are required")
}
if (is.na(opt$ChunkRows) || opt$ChunkRows < 1L) stop("--ChunkRows must be at least 1")
if (!is.finite(opt$ConnectivityZThreshold)) stop("--ConnectivityZThreshold must be finite")
for (path in c(opt$FilteredCalls, opt$RawBed, opt$IntBed, opt$SampleQC)) {
    if (!file.exists(path)) stop("Input file does not exist: ", path)
}

int_bed_paths <- read_file_list(opt$IntBedList, "INT BED")
representative_paths <- read_file_list(opt$RepresentativeList, "Representative CpG")
representatives <- rbindlist(lapply(representative_paths, fread, check.names = FALSE), use.names = TRUE, fill = FALSE)
required_representative_columns <- c("#chr", "start", "end", "phenotype_id", "cluster_id", "cluster_size", "local_connectivity", "selection_type")
if (!identical(names(representatives), required_representative_columns)) {
    stop("Representative-CpG tables must have columns: ", paste(required_representative_columns, collapse = ", "))
}
representatives[, phenotype_id := as.character(phenotype_id)]
if (nrow(representatives) == 0L) stop("No valid CpGs were available as connectivity representatives")
if (anyNA(representatives$phenotype_id) || any(!nzchar(representatives$phenotype_id)) || anyDuplicated(representatives$phenotype_id)) {
    stop("Representative-CpG tables contain missing, blank, or duplicate phenotype IDs")
}
setorder(representatives, `#chr`, start, end, phenotype_id)

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
outliers_output <- paste0(opt$OutputPrefix, ".methylation.connectivity_outliers.tsv")
summary_output <- paste0(opt$OutputPrefix, ".methylation.connectivity_summary.tsv")
representatives_output <- paste0(opt$OutputPrefix, ".methylation.connectivity_representative_cpgs.tsv.gz")
sample_qc_output <- paste0(opt$OutputPrefix, ".methylation.sample_qc.tsv")
filtered_calls_output <- paste0(opt$OutputPrefix, ".methylation.filtered.long.tsv.gz")
raw_bed_output <- paste0(opt$OutputPrefix, ".methylation.raw.bed.gz")
int_bed_output <- paste0(opt$OutputPrefix, ".methylation.INT.bed.gz")
if (any(file.exists(c(outliers_output, summary_output, representatives_output, sample_qc_output, filtered_calls_output, raw_bed_output, int_bed_output)))) {
    stop("One or more connectivity-filter outputs already exist for this prefix")
}

representative_values <- load_representative_values(int_bed_paths, representatives, opt$ChunkRows)
sample_ids <- representative_values$sample_ids
values <- representative_values$values
feature_sd <- apply(values, 1L, sd)
valid_features <- is.finite(feature_sd) & feature_sd > 0
representatives[, used_for_connectivity := valid_features]
fwrite(representatives, representatives_output, sep = "\t", na = "NA")
values <- values[valid_features, , drop = FALSE]

connectivity_score <- rep(NA_real_, length(sample_ids))
connectivity_zscore <- rep(NA_real_, length(sample_ids))
names(connectivity_score) <- sample_ids
names(connectivity_zscore) <- sample_ids
outliers <- data.table(SampleID = character(), Connectivity = numeric(), Z_score = numeric())
method <- "not_computed"

if (nrow(values) < 2L || length(sample_ids) < 3L) {
    message("Not enough nonconstant representative INT CpGs or samples to compute connectivity; keeping all samples")
} else {
    message("Computing full WGCNA sample connectivity from ", nrow(values), " representative INT CpG(s) and ",
            length(sample_ids), " sample(s)")
    normalized_adjacency <- 0.5 + 0.5 * WGCNA::bicor(
        as.data.frame(values, check.names = FALSE),
        use = "pairwise.complete.obs"
    )
    normalized_adjacency[is.na(normalized_adjacency)] <- 0
    connectivity_score <- WGCNA::fundamentalNetworkConcepts(normalized_adjacency)$Connectivity
    names(connectivity_score) <- sample_ids
    connectivity_sd <- sd(connectivity_score, na.rm = TRUE)
    if (!is.finite(connectivity_sd) || connectivity_sd == 0) {
        message("Connectivity scores have zero or undefined variance; keeping all samples")
        method <- "correlation_pruned_full_wgcna_bicor_zero_variance"
    } else {
        connectivity_zscore <- (connectivity_score - mean(connectivity_score, na.rm = TRUE)) / connectivity_sd
        names(connectivity_zscore) <- sample_ids
        outliers <- data.table(
            SampleID = sample_ids[connectivity_zscore < opt$ConnectivityZThreshold],
            Connectivity = connectivity_score[connectivity_zscore < opt$ConnectivityZThreshold],
            Z_score = connectivity_zscore[connectivity_zscore < opt$ConnectivityZThreshold]
        )
        method <- "correlation_pruned_full_wgcna_bicor"
    }
}

retained_samples <- setdiff(sample_ids, outliers$SampleID)
if (length(retained_samples) == 0L) stop("Connectivity filter would remove every sample")
message("Connectivity outlier removal: removed ", nrow(outliers), " of ", length(sample_ids),
        " samples; ", length(retained_samples), " samples remain")

sample_qc <- fread(opt$SampleQC, check.names = FALSE)
if (!("sample_id" %in% names(sample_qc))) stop("Sample-QC file must contain sample_id")
sample_qc[, sample_id := as.character(sample_id)]
if (!setequal(sample_qc$sample_id, sample_ids)) stop("Sample-QC IDs do not match the INT BED sample IDs")
sample_qc[, `:=`(
    connectivity_score = connectivity_score[sample_id],
    connectivity_zscore = connectivity_zscore[sample_id],
    pass_connectivity_filter = sample_id %in% retained_samples
)]
setorder(sample_qc, sample_id)
fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")
fwrite(outliers, outliers_output, sep = "\t", na = "NA")
fwrite(data.table(
    method = method,
    n_representative_cpgs = nrow(representatives),
    n_representative_cpgs_used = nrow(values),
    n_samples_before = length(sample_ids),
    n_samples_removed = nrow(outliers),
    n_samples_retained = length(retained_samples),
    n_samples_in_full_correlation = length(sample_ids),
    z_threshold = opt$ConnectivityZThreshold
), summary_output, sep = "\t")

copy_or_filter_long_calls(opt$FilteredCalls, filtered_calls_output, retained_samples)
copy_or_filter_wide_bed(opt$RawBed, raw_bed_output, retained_samples)
copy_or_filter_wide_bed(opt$IntBed, int_bed_output, retained_samples)

message("Wrote selected representative CpGs: ", representatives_output)
message("Wrote connectivity-filtered long calls: ", filtered_calls_output)
message("Wrote connectivity-filtered raw BED: ", raw_bed_output)
message("Wrote connectivity-filtered INT BED: ", int_bed_output)
