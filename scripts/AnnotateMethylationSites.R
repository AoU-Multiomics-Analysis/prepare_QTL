#!/usr/bin/env Rscript

# Annotate cohort-QC-passing methylation sites with gene/TSS and enhancer context.
# Site and enhancer coordinates are BED (0-based, half-open); GTF coordinates
# are imported in their native 1-based closed representation.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(plyranges)
})

option_list <- list(
    make_option("--SiteMetadata", type = "character", help = "All-site methylation metadata TSV(.gz) [required]"),
    make_option("--AnnotationGTF", type = "character", help = "Gene annotation GTF(.gz) [required]"),
    make_option("--EnhancerAnnotations", type = "character", help = "Six-column enhancer BED-like file [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]"),
    make_option("--PromoterWindow", type = "integer", default = 2000,
                help = "Bases on either side of a TSS defining a promoter [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("SiteMetadata", "AnnotationGTF", "EnhancerAnnotations", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--SiteMetadata, --AnnotationGTF, --EnhancerAnnotations, and --OutputPrefix are required")
}
for (path in unlist(opt[required_options[required_options != "OutputPrefix"]])) {
    if (!file.exists(path)) stop("Input file does not exist: ", path)
}
if (is.na(opt$PromoterWindow) || opt$PromoterWindow < 0) stop("--PromoterWindow must be a non-negative integer")

collapse_values <- function(x) {
    values <- unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))
    if (length(values) == 0) NA_character_ else paste(sort(values), collapse = ";")
}

summarize_overlaps <- function(query_gr, subject_gr, fields) {
    hits <- GenomicRanges::findOverlaps(query_gr, subject_gr, ignore.strand = TRUE)
    if (length(hits) == 0) {
        result <- data.table(site_index = integer(), n_overlaps = integer())
        for (field in fields) result[, (field) := character()]
        return(result)
    }
    result <- data.table(site_index = S4Vectors::queryHits(hits))
    metadata <- as.data.table(as.data.frame(S4Vectors::mcols(subject_gr)))
    result[, (fields) := metadata[S4Vectors::subjectHits(hits), fields, with = FALSE]]
    result[, c(list(n_overlaps = .N), lapply(.SD, collapse_values)), by = site_index, .SDcols = fields]
}

add_overlap_columns <- function(annotation, overlaps, count_column, fields) {
    annotation[, (count_column) := 0L]
    for (field in fields) annotation[, (field) := NA_character_]
    if (nrow(overlaps) == 0) return(annotation)
    matched <- match(annotation$site_index, overlaps$site_index)
    found <- !is.na(matched)
    annotation[found, (count_column) := as.integer(overlaps$n_overlaps[matched[found]])]
    for (field in fields) annotation[found, (field) := overlaps[[field]][matched[found]]]
    annotation
}

metadata <- fread(opt$SiteMetadata)
required_site_columns <- c("#chrom", "begin", "end", "site_key", "keep_site")
missing_site_columns <- setdiff(required_site_columns, names(metadata))
if (length(missing_site_columns) > 0) stop("Site metadata is missing: ", paste(missing_site_columns, collapse = ", "))
sites <- metadata[keep_site == TRUE]
if (anyNA(sites$begin) || anyNA(sites$end) || any(sites$begin < 0) || any(sites$end <= sites$begin)) {
    stop("Passing site metadata must contain valid BED start/end coordinates")
}
output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- paste0(opt$OutputPrefix, ".methylation.passing_site_annotations.tsv.gz")
message("Annotating ", nrow(sites), " cohort-QC-passing methylation sites")
if (nrow(sites) == 0) {
    fwrite(sites, output_file, sep = "\t", na = "NA")
    message("No passing sites to annotate; wrote empty annotation file: ", output_file)
    quit(save = "no", status = 0)
}

# as_granges() from plyranges constructs GRanges while retaining metadata.
sites[, site_index := seq_len(.N)]
site_gr <- plyranges::as_granges(data.frame(
    seqnames = as.character(sites[["#chrom"]]), start = as.integer(sites$begin) + 1L,
    end = as.integer(sites$end), strand = "*", site_index = sites$site_index
))

message("Loading GTF gene annotations")
gtf <- rtracklayer::import(opt$AnnotationGTF)
gtf_metadata <- S4Vectors::mcols(gtf)
if (!("type" %in% names(gtf_metadata))) stop("GTF does not contain a feature 'type' column")
genes <- gtf[as.character(gtf_metadata$type) == "gene"]
if (length(genes) == 0) stop("GTF contains no features with type == 'gene'")
gene_metadata <- S4Vectors::mcols(genes)
if (!("gene_id" %in% names(gene_metadata))) stop("GTF gene features must contain a gene_id attribute")
gene_ids <- as.character(gene_metadata$gene_id)
gene_names <- if ("gene_name" %in% names(gene_metadata)) as.character(gene_metadata$gene_name) else rep(NA_character_, length(genes))
gene_strands <- as.character(GenomicRanges::strand(genes))
gene_table <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(genes)), start = GenomicRanges::start(genes),
    end = GenomicRanges::end(genes), strand = gene_strands, gene_id = gene_ids, gene_name = gene_names
)
gene_gr <- plyranges::as_granges(gene_table)
tss_position <- ifelse(gene_strands == "-", GenomicRanges::end(genes), GenomicRanges::start(genes))
tss_table <- data.frame(
    seqnames = gene_table$seqnames, start = tss_position, end = tss_position,
    strand = gene_strands, gene_id = gene_ids, gene_name = gene_names
)
tss_gr <- plyranges::as_granges(tss_table)
promoter_table <- tss_table
promoter_table$start <- pmax(1L, promoter_table$start - opt$PromoterWindow)
promoter_table$end <- promoter_table$end + opt$PromoterWindow
promoter_gr <- plyranges::as_granges(promoter_table)

nearest_tss <- GenomicRanges::nearest(site_gr, tss_gr, ignore.strand = TRUE)
sites[, `:=`(
    nearest_tss_distance = NA_integer_, nearest_tss_gene_id = NA_character_,
    nearest_tss_gene_name = NA_character_, nearest_tss_strand = NA_character_
)]
has_nearest <- !is.na(nearest_tss)
if (any(has_nearest)) {
    queries <- which(has_nearest)
    subjects <- nearest_tss[has_nearest]
    sites[queries, `:=`(
        nearest_tss_distance = as.integer(abs(GenomicRanges::start(site_gr)[queries] - GenomicRanges::start(tss_gr)[subjects])),
        nearest_tss_gene_id = S4Vectors::mcols(tss_gr)$gene_id[subjects],
        nearest_tss_gene_name = S4Vectors::mcols(tss_gr)$gene_name[subjects],
        nearest_tss_strand = as.character(GenomicRanges::strand(tss_gr)[subjects])
    )]
}

message("Annotating promoter and gene-body overlaps")
promoter_hits <- summarize_overlaps(site_gr, promoter_gr, c("gene_id", "gene_name"))
sites <- add_overlap_columns(sites, promoter_hits, "n_promoter_genes", c("gene_id", "gene_name"))
setnames(sites, c("gene_id", "gene_name"), c("promoter_gene_id", "promoter_gene_name"))
sites[, in_promoter := n_promoter_genes > 0]
gene_body_hits <- summarize_overlaps(site_gr, gene_gr, c("gene_id", "gene_name"))
sites <- add_overlap_columns(sites, gene_body_hits, "n_gene_body_genes", c("gene_id", "gene_name"))
setnames(sites, c("gene_id", "gene_name"), c("gene_body_gene_id", "gene_body_gene_name"))
sites[, in_gene_body := n_gene_body_genes > 0]
sites[, genomic_context := fcase(in_promoter, "promoter", in_gene_body, "gene_body", default = "intergenic")]

message("Loading and annotating enhancer overlaps")
enhancers <- fread(opt$EnhancerAnnotations, header = FALSE)
if (ncol(enhancers) < 6) stop("Enhancer annotations must have at least six columns (V1-V6)")
setnames(enhancers, names(enhancers)[seq_len(6)], c("enhancer_chrom", "enhancer_begin", "enhancer_end", "enhancer_v4_id", "enhancer_v5_id", "enhancer_type"))
enhancers[, `:=`(enhancer_begin = as.integer(enhancer_begin), enhancer_end = as.integer(enhancer_end))]
if (anyNA(enhancers$enhancer_begin) || anyNA(enhancers$enhancer_end) || any(enhancers$enhancer_begin < 0) || any(enhancers$enhancer_end <= enhancers$enhancer_begin)) {
    stop("Enhancer annotations must use valid BED start/end coordinates")
}
enhancer_gr <- plyranges::as_granges(data.frame(
    seqnames = as.character(enhancers$enhancer_chrom), start = enhancers$enhancer_begin + 1L,
    end = enhancers$enhancer_end, strand = "*", enhancer_v4_id = as.character(enhancers$enhancer_v4_id),
    enhancer_v5_id = as.character(enhancers$enhancer_v5_id), enhancer_type = as.character(enhancers$enhancer_type)
))
enhancer_hits <- summarize_overlaps(site_gr, enhancer_gr, c("enhancer_v4_id", "enhancer_v5_id", "enhancer_type"))
sites <- add_overlap_columns(sites, enhancer_hits, "n_overlapping_enhancers", c("enhancer_v4_id", "enhancer_v5_id", "enhancer_type"))
sites[, in_enhancer := n_overlapping_enhancers > 0]

setorder(sites, `#chrom`, begin, end)
sites[, site_index := NULL]
fwrite(sites, output_file, sep = "\t", na = "NA")
message("Wrote annotations for ", nrow(sites), " passing sites: ", output_file)
