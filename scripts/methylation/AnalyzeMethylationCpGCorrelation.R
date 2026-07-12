#!/usr/bin/env Rscript

# Stream a methylation phenotype BED to characterize local CpG correlation
# clusters without loading a whole chromosome-wide phenotype matrix in memory.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(optparse)
})

is_present <- function(x) {
    !is.null(x) && !is.na(x) && nzchar(x) && !tolower(x) %in% c("none", "null", "na")
}

make_empty_clusters <- function() {
    data.table(
        cluster_id = integer(), `#chr` = character(), cluster_start = integer(), cluster_end = integer(),
        cluster_span_bp = integer(), n_correlated_cpgs = integer(), n_correlated_pairs = integer(),
        mean_pair_distance_bp = numeric(), max_pair_distance_bp = integer(),
        mean_abs_correlation = numeric(), max_abs_correlation = numeric()
    )
}

load_covariate_design <- function(covariate_path, bed_sample_ids) {
    covariate_table <- fread(covariate_path, header = FALSE, check.names = FALSE,
                             na.strings = c("", "NA", "NaN", "."))
    if (nrow(covariate_table) < 2 || ncol(covariate_table) < 2) {
        stop("--Covariates must use TensorQTL orientation: one covariate per row and one sample per column")
    }

    covariate_sample_ids <- as.character(unlist(covariate_table[1, -1, with = FALSE]))
    covariate_names <- as.character(covariate_table[[1]][-1])
    if (anyNA(covariate_sample_ids) || any(!nzchar(covariate_sample_ids)) || anyDuplicated(covariate_sample_ids)) {
        stop("Covariates header contains missing, blank, or duplicated sample IDs")
    }
    if (anyNA(covariate_names) || any(!nzchar(covariate_names))) stop("Covariates contain a missing or blank covariate ID")

    raw_values <- as.matrix(covariate_table[-1, -1, with = FALSE])
    numeric_values <- suppressWarnings(matrix(as.numeric(raw_values), nrow = nrow(raw_values)))
    if (any(is.na(numeric_values) & !is.na(raw_values))) stop("Covariates contain non-numeric values")
    covariate_by_sample <- as.data.frame(t(numeric_values), stringsAsFactors = FALSE)
    rownames(covariate_by_sample) <- covariate_sample_ids
    colnames(covariate_by_sample) <- make.names(make.unique(covariate_names), unique = TRUE)

    shared_samples <- bed_sample_ids[bed_sample_ids %in% rownames(covariate_by_sample)]
    if (length(shared_samples) < 3) stop("Fewer than three BED samples overlap the covariates")
    if (length(shared_samples) < length(bed_sample_ids)) {
        message(length(bed_sample_ids) - length(shared_samples), " BED sample(s) are absent from covariates and will be excluded")
    }
    covariate_by_sample <- covariate_by_sample[shared_samples, , drop = FALSE]
    complete_samples <- complete.cases(covariate_by_sample)
    if (!all(complete_samples)) {
        message(sum(!complete_samples), " sample(s) are excluded because covariates are missing")
        covariate_by_sample <- covariate_by_sample[complete_samples, , drop = FALSE]
        shared_samples <- rownames(covariate_by_sample)
    }
    design <- model.matrix(~ ., data = covariate_by_sample)
    design_qr <- qr(design)
    if (nrow(design) <= design_qr$rank) stop("Not enough samples after covariate filtering to residualize CpGs")

    list(
        sample_ids = shared_samples,
        qr = design_qr,
        n_covariates = ncol(design) - 1L
    )
}

residualize_and_standardize <- function(values, design_qr = NULL) {
    if (any(!is.finite(values))) stop("Input BED contains missing or non-finite methylation values; use the imputed cohort BED")
    if (!is.null(design_qr)) {
        values <- t(qr.resid(design_qr, t(values)))
    }
    centers <- rowMeans(values)
    centered <- values - centers
    norms <- sqrt(rowSums(centered^2))
    valid <- is.finite(norms) & norms > 0
    standardized <- matrix(NA_real_, nrow = nrow(values), ncol = ncol(values))
    if (any(valid)) standardized[valid, ] <- centered[valid, , drop = FALSE] / norms[valid]
    list(values = standardized, valid = valid)
}

open_bed_connection <- function(path) {
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        pipe(sprintf("bgzip -c -d %s", shQuote(path)), open = "r")
    } else {
        file(path, open = "r")
    }
}

plot_empty <- function(path, title, subtitle) {
    plot <- ggplot() +
        annotate("text", x = 0, y = 0, label = subtitle, size = 5) +
        xlim(-1, 1) + ylim(-1, 1) +
        labs(title = title) +
        theme_void()
    ggsave(path, plot, width = 8, height = 5, dpi = 200, bg = "white")
}

option_list <- list(
    make_option("--InputBed", type = "character", help = "Cohort INT methylation BED (optionally bgzipped) [required]"),
    make_option("--OutputPrefix", type = "character", help = "Output prefix [required]"),
    make_option("--Covariates", type = "character", default = NULL,
                help = "Optional TensorQTL-format covariates TSV; CpGs are residualized before correlation"),
    make_option("--WindowBP", type = "integer", default = 1000,
                help = "Maximum distance between CpGs considered a local pair [default: %default]"),
    make_option("--MinAbsCorrelation", type = "double", default = 0.95,
                help = "Absolute Pearson-correlation threshold for a correlated pair [default: %default]"),
    make_option("--ChunkRows", type = "integer", default = 1000,
                help = "Number of CpG rows parsed per streaming block [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!is_present(opt$InputBed) || !is_present(opt$OutputPrefix)) stop("--InputBed and --OutputPrefix are required")
if (!file.exists(opt$InputBed)) stop("Input BED does not exist: ", opt$InputBed)
if (!is.null(opt$Covariates) && !file.exists(opt$Covariates)) stop("Covariates file does not exist: ", opt$Covariates)
if (is.na(opt$WindowBP) || opt$WindowBP < 1) stop("--WindowBP must be at least 1")
if (!is.finite(opt$MinAbsCorrelation) || opt$MinAbsCorrelation <= 0 || opt$MinAbsCorrelation > 1) {
    stop("--MinAbsCorrelation must be in (0, 1]")
}
if (is.na(opt$ChunkRows) || opt$ChunkRows < 1) stop("--ChunkRows must be at least 1")

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cluster_table_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_clusters.tsv")
cluster_output_path <- paste0(cluster_table_path, ".gz")
summary_output_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_summary.tsv")
cluster_size_plot_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_cluster_sizes.png")
mean_distance_plot_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_mean_distances.png")
max_distance_plot_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_max_distances.png")
span_plot_path <- paste0(opt$OutputPrefix, ".methylation.cpg_correlation_span_vs_size.png")
if (file.exists(cluster_table_path) || file.exists(cluster_output_path)) {
    stop("Correlation-cluster output already exists for this prefix: ", cluster_output_path)
}

con <- open_bed_connection(opt$InputBed)
on.exit(close(con), add = TRUE)
header_line <- readLines(con, n = 1L, warn = FALSE)
if (length(header_line) != 1L) stop("Input BED is empty")
header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
if (length(header) < 5L) stop("Input BED must contain four BED metadata columns and at least one sample column")
metadata_columns <- header[seq_len(4L)]
bed_sample_ids <- header[-seq_len(4L)]
if (anyNA(bed_sample_ids) || any(!nzchar(bed_sample_ids)) || anyDuplicated(bed_sample_ids)) {
    stop("Input BED contains missing, blank, or duplicated sample IDs")
}

if (is_present(opt$Covariates)) {
    covariate_design <- load_covariate_design(opt$Covariates, bed_sample_ids)
    sample_ids <- covariate_design$sample_ids
    sample_indices <- match(sample_ids, bed_sample_ids)
    message("Residualizing CpGs with ", covariate_design$n_covariates, " covariate(s) across ", length(sample_ids), " sample(s)")
} else {
    covariate_design <- NULL
    sample_ids <- bed_sample_ids
    sample_indices <- seq_along(bed_sample_ids)
    message("No covariates supplied; using Pearson correlation of centered/scaled methylation values across ", length(sample_ids), " sample(s)")
}

state <- new.env(parent = emptyenv())
state$current_chromosome <- NULL
state$last_start <- -Inf
state$active_start <- numeric()
state$active_site <- character()
state$active_component <- integer()
state$active_values <- list()
state$active_head <- 1L
state$components <- new.env(parent = emptyenv())
state$next_component_id <- 1L
state$next_cluster_id <- 1L
state$cluster_buffer <- list()
state$clusters_written <- FALSE
state$n_cpgs <- 0L
state$n_zero_variance_cpgs <- 0L
state$n_correlated_pairs <- 0L

component_key <- function(component_id) as.character(component_id)

get_component <- function(component_id) state$components[[component_key(component_id)]]
set_component <- function(component_id, component) state$components[[component_key(component_id)]] <- component

new_component <- function(chromosome) {
    component_id <- state$next_component_id
    state$next_component_id <- state$next_component_id + 1L
    set_component(component_id, list(
        chromosome = chromosome, start = Inf, end = -Inf, n_cpgs = 0L, active_members = 0L,
        n_pairs = 0L, sum_distance = 0, max_distance = 0, sum_abs_correlation = 0, max_abs_correlation = 0
    ))
    component_id
}

flush_cluster_buffer <- function() {
    if (length(state$cluster_buffer) == 0L) return(invisible(NULL))
    cluster_table <- rbindlist(state$cluster_buffer, use.names = TRUE)
    fwrite(cluster_table, cluster_table_path, sep = "\t", append = state$clusters_written)
    state$clusters_written <- TRUE
    state$cluster_buffer <- list()
    invisible(NULL)
}

finalize_component <- function(component_id) {
    component <- get_component(component_id)
    if (is.null(component)) return(invisible(NULL))
    if (component$n_pairs > 0L) {
        state$cluster_buffer[[length(state$cluster_buffer) + 1L]] <- data.table(
            cluster_id = state$next_cluster_id,
            `#chr` = component$chromosome,
            cluster_start = as.integer(component$start),
            cluster_end = as.integer(component$end),
            cluster_span_bp = as.integer(component$end - component$start),
            n_correlated_cpgs = as.integer(component$n_cpgs),
            n_correlated_pairs = as.integer(component$n_pairs),
            mean_pair_distance_bp = component$sum_distance / component$n_pairs,
            max_pair_distance_bp = as.integer(component$max_distance),
            mean_abs_correlation = component$sum_abs_correlation / component$n_pairs,
            max_abs_correlation = component$max_abs_correlation
        )
        state$next_cluster_id <- state$next_cluster_id + 1L
        if (length(state$cluster_buffer) >= 1000L) flush_cluster_buffer()
    }
    rm(list = component_key(component_id), envir = state$components)
    invisible(NULL)
}

merge_components <- function(target_id, source_id) {
    if (target_id == source_id) return(target_id)
    target <- get_component(target_id)
    source <- get_component(source_id)
    if (is.null(source)) return(target_id)
    target$start <- min(target$start, source$start)
    target$end <- max(target$end, source$end)
    target$n_cpgs <- target$n_cpgs + source$n_cpgs
    target$active_members <- target$active_members + source$active_members
    target$n_pairs <- target$n_pairs + source$n_pairs
    target$sum_distance <- target$sum_distance + source$sum_distance
    target$max_distance <- max(target$max_distance, source$max_distance)
    target$sum_abs_correlation <- target$sum_abs_correlation + source$sum_abs_correlation
    target$max_abs_correlation <- max(target$max_abs_correlation, source$max_abs_correlation)
    set_component(target_id, target)
    state$active_component[state$active_component == source_id] <- target_id
    rm(list = component_key(source_id), envir = state$components)
    target_id
}

add_existing_active_node <- function(component_id, active_index) {
    component <- get_component(component_id)
    component$n_cpgs <- component$n_cpgs + 1L
    component$active_members <- component$active_members + 1L
    component$start <- min(component$start, state$active_start[active_index])
    component$end <- max(component$end, state$active_start[active_index])
    set_component(component_id, component)
    state$active_component[active_index] <- component_id
}

add_new_node <- function(component_id, chromosome, start) {
    component <- get_component(component_id)
    component$n_cpgs <- component$n_cpgs + 1L
    component$active_members <- component$active_members + 1L
    component$start <- min(component$start, start)
    component$end <- max(component$end, start)
    set_component(component_id, component)
}

add_edges <- function(component_id, distances, correlations) {
    component <- get_component(component_id)
    component$n_pairs <- component$n_pairs + length(distances)
    component$sum_distance <- component$sum_distance + sum(distances)
    component$max_distance <- max(component$max_distance, max(distances))
    component$sum_abs_correlation <- component$sum_abs_correlation + sum(abs(correlations))
    component$max_abs_correlation <- max(component$max_abs_correlation, max(abs(correlations)))
    set_component(component_id, component)
}

expire_active_until <- function(minimum_start) {
    active_length <- length(state$active_start)
    while (state$active_head <= active_length && state$active_start[state$active_head] < minimum_start) {
        component_id <- state$active_component[state$active_head]
        if (component_id > 0L) {
            component <- get_component(component_id)
            component$active_members <- component$active_members - 1L
            set_component(component_id, component)
            if (component$active_members == 0L) finalize_component(component_id)
        }
        state$active_head <- state$active_head + 1L
    }
    if (state$active_head > active_length) {
        state$active_start <- numeric()
        state$active_site <- character()
        state$active_component <- integer()
        state$active_values <- list()
        state$active_head <- 1L
    } else if (state$active_head > 1024L) {
        retained <- seq.int(state$active_head, active_length)
        state$active_start <- state$active_start[retained]
        state$active_site <- state$active_site[retained]
        state$active_component <- state$active_component[retained]
        state$active_values <- state$active_values[retained]
        state$active_head <- 1L
    }
}

flush_chromosome <- function() {
    if (is.null(state$current_chromosome)) return(invisible(NULL))
    active_length <- length(state$active_start)
    if (state$active_head <= active_length) {
        for (active_index in seq.int(state$active_head, active_length)) {
            component_id <- state$active_component[active_index]
            if (component_id > 0L) {
                component <- get_component(component_id)
                component$active_members <- component$active_members - 1L
                set_component(component_id, component)
                if (component$active_members == 0L) finalize_component(component_id)
            }
        }
    }
    state$active_start <- numeric()
    state$active_site <- character()
    state$active_component <- integer()
    state$active_values <- list()
    state$active_head <- 1L
    state$current_chromosome <- NULL
    state$last_start <- -Inf
    invisible(NULL)
}

process_cpg <- function(chromosome, start, site_id, standardized_values, is_valid) {
    if (is.null(state$current_chromosome) || chromosome != state$current_chromosome) {
        flush_chromosome()
        state$current_chromosome <- chromosome
    }
    if (start < state$last_start) stop("Input BED must be sorted by chromosome and start coordinate")
    state$last_start <- start
    state$n_cpgs <- state$n_cpgs + 1L
    if (!is_valid) {
        state$n_zero_variance_cpgs <- state$n_zero_variance_cpgs + 1L
        return(invisible(NULL))
    }

    expire_active_until(start - opt$WindowBP)
    active_length <- length(state$active_start)
    candidate_indices <- if (state$active_head <= active_length) seq.int(state$active_head, active_length) else integer()
    matching_indices <- integer()
    matching_correlations <- numeric()
    if (length(candidate_indices) > 0L) {
        candidate_values <- do.call(rbind, state$active_values[candidate_indices])
        correlations <- drop(candidate_values %*% standardized_values)
        matches <- abs(correlations) >= opt$MinAbsCorrelation
        matching_indices <- candidate_indices[matches]
        matching_correlations <- correlations[matches]
    }

    component_id <- 0L
    if (length(matching_indices) > 0L) {
        existing_components <- unique(state$active_component[matching_indices])
        existing_components <- existing_components[existing_components > 0L]
        component_id <- if (length(existing_components) == 0L) new_component(chromosome) else existing_components[[1]]
        if (length(existing_components) > 1L) {
            for (source_id in existing_components[-1L]) component_id <- merge_components(component_id, source_id)
        }
        for (active_index in matching_indices) {
            active_component_id <- state$active_component[active_index]
            if (active_component_id == 0L) {
                add_existing_active_node(component_id, active_index)
            } else if (active_component_id != component_id) {
                component_id <- merge_components(component_id, active_component_id)
            }
        }
        add_new_node(component_id, chromosome, start)
        distances <- start - state$active_start[matching_indices]
        add_edges(component_id, distances, matching_correlations)
        state$n_correlated_pairs <- state$n_correlated_pairs + length(matching_indices)
    }

    state$active_start <- c(state$active_start, start)
    state$active_site <- c(state$active_site, site_id)
    state$active_component <- c(state$active_component, component_id)
    state$active_values[[length(state$active_values) + 1L]] <- standardized_values
    invisible(NULL)
}

message("Streaming BED: ", opt$InputBed)
repeat {
    lines <- readLines(con, n = opt$ChunkRows, warn = FALSE)
    if (length(lines) == 0L) break
    block <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), header = FALSE, sep = "\t", col.names = header,
                   check.names = FALSE, showProgress = FALSE)
    if (ncol(block) != length(header)) stop("BED row has a different number of columns than the header")
    starts <- suppressWarnings(as.integer(block[[metadata_columns[[2]]]]))
    if (anyNA(starts)) stop("BED start column is not integer-like")
    chromosomes <- as.character(block[[metadata_columns[[1]]]])
    site_ids <- as.character(block[[metadata_columns[[4]]]])
    raw_values <- as.matrix(block[, ..bed_sample_ids])
    values <- suppressWarnings(matrix(as.numeric(raw_values), nrow = nrow(block), dimnames = list(NULL, bed_sample_ids)))
    if (any(is.na(values) & !is.na(raw_values))) stop("BED phenotype values are not numeric")
    values <- values[, sample_indices, drop = FALSE]
    standardized <- residualize_and_standardize(values, if (is.null(covariate_design)) NULL else covariate_design$qr)
    for (row_index in seq_len(nrow(block))) {
        process_cpg(chromosomes[[row_index]], starts[[row_index]], site_ids[[row_index]],
                    standardized$values[row_index, ], standardized$valid[[row_index]])
    }
}
flush_chromosome()
flush_cluster_buffer()
if (!state$clusters_written) fwrite(make_empty_clusters(), cluster_table_path, sep = "\t")
status <- if (nzchar(Sys.which("bgzip"))) system2("bgzip", c("-f", cluster_table_path)) else 1L
if (!identical(status, 0L)) {
    if (nzchar(Sys.which("gzip"))) {
        status <- system2("gzip", c("-f", cluster_table_path))
    }
    if (!identical(status, 0L)) stop("Failed to compress correlation-cluster output")
}

cluster_table <- fread(cluster_output_path)
summary_table <- data.table(
    input_bed = opt$InputBed,
    covariates = if (is_present(opt$Covariates)) opt$Covariates else "",
    n_samples_used = length(sample_ids),
    n_cpgs_processed = state$n_cpgs,
    n_zero_variance_cpgs = state$n_zero_variance_cpgs,
    window_bp = opt$WindowBP,
    min_abs_correlation = opt$MinAbsCorrelation,
    n_correlated_pairs = state$n_correlated_pairs,
    n_clusters = nrow(cluster_table),
    n_cpgs_in_clusters = if (nrow(cluster_table) == 0L) 0L else sum(cluster_table$n_correlated_cpgs),
    mean_cluster_size = if (nrow(cluster_table) == 0L) NA_real_ else mean(cluster_table$n_correlated_cpgs),
    median_cluster_size = if (nrow(cluster_table) == 0L) NA_real_ else median(cluster_table$n_correlated_cpgs),
    mean_correlated_pair_distance_bp = if (nrow(cluster_table) == 0L) NA_real_ else weighted.mean(cluster_table$mean_pair_distance_bp, cluster_table$n_correlated_pairs),
    max_correlated_pair_distance_bp = if (nrow(cluster_table) == 0L) NA_integer_ else max(cluster_table$max_pair_distance_bp)
)
fwrite(summary_table, summary_output_path, sep = "\t")

if (nrow(cluster_table) == 0L) {
    plot_empty(cluster_size_plot_path, "Correlated CpG cluster sizes", "No local CpG pairs met the correlation threshold")
    plot_empty(mean_distance_plot_path, "Mean correlated-pair distance", "No local CpG pairs met the correlation threshold")
    plot_empty(max_distance_plot_path, "Maximum correlated-pair distance", "No local CpG pairs met the correlation threshold")
    plot_empty(span_plot_path, "CpG cluster span versus size", "No local CpG pairs met the correlation threshold")
} else {
    cluster_size_plot <- ggplot(cluster_table, aes(x = n_correlated_cpgs)) +
        geom_histogram(binwidth = 1, boundary = 0.5, fill = "#2C7FB8", color = "white") +
        labs(title = "Local correlated CpG cluster sizes", x = "Correlated CpGs per cluster", y = "Number of clusters") +
        theme_minimal(base_size = 12)
    mean_distance_plot <- ggplot(cluster_table, aes(x = mean_pair_distance_bp)) +
        geom_histogram(bins = 50, fill = "#41AB5D", color = "white") +
        labs(title = "Mean distance between correlated CpGs", x = "Mean correlated-pair distance (bp)", y = "Number of clusters") +
        theme_minimal(base_size = 12)
    max_distance_plot <- ggplot(cluster_table, aes(x = max_pair_distance_bp)) +
        geom_histogram(bins = 50, fill = "#8856A7", color = "white") +
        labs(title = "Maximum direct-pair distance in each cluster", x = "Maximum correlated-pair distance (bp)", y = "Number of clusters") +
        theme_minimal(base_size = 12)
    span_plot <- ggplot(cluster_table, aes(x = cluster_span_bp, y = n_correlated_cpgs)) +
        geom_point(alpha = 0.35, color = "#D95F0E") +
        scale_x_continuous(trans = "log1p") +
        scale_y_continuous(trans = "log1p") +
        labs(title = "Local correlated CpG cluster span and size", x = "Cluster span (bp; log1p scale)", y = "Correlated CpGs (log1p scale)") +
        theme_minimal(base_size = 12)
    ggsave(cluster_size_plot_path, cluster_size_plot, width = 8, height = 5, dpi = 200, bg = "white")
    ggsave(mean_distance_plot_path, mean_distance_plot, width = 8, height = 5, dpi = 200, bg = "white")
    ggsave(max_distance_plot_path, max_distance_plot, width = 8, height = 5, dpi = 200, bg = "white")
    ggsave(span_plot_path, span_plot, width = 8, height = 5, dpi = 200, bg = "white")
}

message("Processed ", state$n_cpgs, " CpGs across ", length(sample_ids), " sample(s)")
message("Detected ", state$n_correlated_pairs, " correlated local pairs in ", nrow(cluster_table), " cluster(s)")
message("Wrote cluster table: ", cluster_output_path)
message("Wrote summary: ", summary_output_path)
