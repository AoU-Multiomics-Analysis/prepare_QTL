#!/usr/bin/env Rscript

# Annotate cohort-QC-passing methylation sites with gene/TSS, cCRE, and CpG-island context.
# Site, cCRE, and CpG-island coordinates are BED (0-based, half-open); GTF coordinates
# are imported in their native 1-based closed representation.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(plyranges)
})

option_list <- list(
    make_option("--SiteMetadata", type = "character", help = "All-site methylation metadata TSV(.gz) [required]"),
    make_option("--AnnotationGTF", type = "character", help = "Gene annotation GTF(.gz) [required]"),
    make_option("--CCREAnnotations", type = "character", help = "Six-column ENCODE cCRE BED-like file [required]"),
    make_option("--CpGIslandAnnotations", type = "character", help = "UCSC cpgIslandExt table [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]"),
    make_option("--PromoterWindow", type = "integer", default = 2000,
                help = "Bases on either side of a TSS defining a promoter [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("SiteMetadata", "AnnotationGTF", "CCREAnnotations", "CpGIslandAnnotations", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--SiteMetadata, --AnnotationGTF, --CCREAnnotations, --CpGIslandAnnotations, and --OutputPrefix are required")
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
    if (is.null(subject_gr)) {
        result <- data.table(site_index = integer(), n_overlaps = integer())
        for (field in fields) result[, (field) := character()]
        return(result)
    }
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

annotate_gene_feature <- function(annotation, query_gr, feature_gr, label) {
    fields <- c("gene_id", "gene_name")
    overlaps <- summarize_overlaps(query_gr, feature_gr, fields)
    annotation <- add_overlap_columns(annotation, overlaps, paste0("n_", label, "_genes"), fields)
    setnames(annotation, fields, paste0(label, c("_gene_id", "_gene_name")))
    annotation[, (paste0("in_", label)) := get(paste0("n_", label, "_genes")) > 0]
    annotation
}

split_collapsed_values <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1]]
}

derive_intron_overlaps <- function(gene_body_hits, exon_hits, gene_lookup) {
    if (nrow(gene_body_hits) == 0) {
        return(data.table(site_index = integer(), n_overlaps = integer(), gene_id = character(), gene_name = character()))
    }
    exon_by_site <- setNames(exon_hits$gene_id, exon_hits$site_index)
    intron_rows <- lapply(seq_len(nrow(gene_body_hits)), function(i) {
        body_ids <- split_collapsed_values(gene_body_hits$gene_id[i])
        exon_ids <- split_collapsed_values(exon_by_site[as.character(gene_body_hits$site_index[i])])
        intron_ids <- setdiff(body_ids, exon_ids)
        if (length(intron_ids) == 0) return(NULL)
        intron_names <- gene_lookup[match(intron_ids, gene_id), gene_name]
        data.table(
            site_index = gene_body_hits$site_index[i],
            n_overlaps = length(intron_ids),
            gene_id = collapse_values(intron_ids),
            gene_name = collapse_values(intron_names)
        )
    })
    intron_rows <- Filter(Negate(is.null), intron_rows)
    if (length(intron_rows) == 0) {
        return(data.table(site_index = integer(), n_overlaps = integer(), gene_id = character(), gene_name = character()))
    }
    rbindlist(intron_rows, use.names = TRUE)
}

read_cpg_islands <- function(file_path) {
    islands <- fread(file_path)
    required_columns <- c("chrom", "chromStart", "chromEnd", "name")
    if (!all(required_columns %in% names(islands))) {
        islands <- fread(file_path, header = FALSE)
        if (ncol(islands) < 5) stop("UCSC cpgIslandExt input must have at least five columns")
        schema <- c("bin", "chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
        setnames(islands, names(islands)[seq_len(min(ncol(islands), length(schema)))], schema[seq_len(min(ncol(islands), length(schema)))])
    }
    missing_columns <- setdiff(required_columns, names(islands))
    if (length(missing_columns) > 0) stop("CpG-island input is missing: ", paste(missing_columns, collapse = ", "))
    islands[, `:=`(chromStart = as.integer(chromStart), chromEnd = as.integer(chromEnd))]
    if (anyNA(islands$chromStart) || anyNA(islands$chromEnd) || any(islands$chromStart < 0) || any(islands$chromEnd <= islands$chromStart)) {
        stop("CpG-island input must use valid UCSC BED chromStart/chromEnd coordinates")
    }
    islands
}

expand_ranges <- function(granges, bases) {
    plyranges::as_granges(data.frame(
        seqnames = as.character(GenomicRanges::seqnames(granges)),
        start = pmax(1L, GenomicRanges::start(granges) - bases),
        end = GenomicRanges::end(granges) + bases,
        strand = "*"
    ))
}

overlaps_any <- function(query_gr, subject_gr) {
    hits <- GenomicRanges::findOverlaps(query_gr, subject_gr, ignore.strand = TRUE)
    present <- rep(FALSE, length(query_gr))
    present[unique(S4Vectors::queryHits(hits))] <- TRUE
    present
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
make_feature_gr <- function(feature_types) {
    features <- gtf[as.character(gtf_metadata$type) %chin% feature_types]
    if (length(features) == 0) return(NULL)
    feature_metadata <- S4Vectors::mcols(features)
    if (!("gene_id" %in% names(feature_metadata))) {
        stop("GTF ", paste(feature_types, collapse = "/"), " features must contain a gene_id attribute")
    }
    feature_names <- if ("gene_name" %in% names(feature_metadata)) as.character(feature_metadata$gene_name) else rep(NA_character_, length(features))
    plyranges::as_granges(data.frame(
        seqnames = as.character(GenomicRanges::seqnames(features)),
        start = GenomicRanges::start(features),
        end = GenomicRanges::end(features),
        strand = as.character(GenomicRanges::strand(features)),
        gene_id = as.character(feature_metadata$gene_id),
        gene_name = feature_names
    ))
}
exon_gr <- make_feature_gr("exon")
cds_gr <- make_feature_gr("CDS")
five_prime_utr_gr <- make_feature_gr(c("five_prime_utr", "5UTR"))
three_prime_utr_gr <- make_feature_gr(c("three_prime_utr", "3UTR"))
utr_gr <- make_feature_gr("UTR")
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

message("Annotating promoter, gene-body, and GTF subfeature overlaps")
sites <- annotate_gene_feature(sites, site_gr, promoter_gr, "promoter")
gene_body_hits <- summarize_overlaps(site_gr, gene_gr, c("gene_id", "gene_name"))
sites <- add_overlap_columns(sites, gene_body_hits, "n_gene_body_genes", c("gene_id", "gene_name"))
setnames(sites, c("gene_id", "gene_name"), c("gene_body_gene_id", "gene_body_gene_name"))
sites[, in_gene_body := n_gene_body_genes > 0]
exon_hits <- summarize_overlaps(site_gr, exon_gr, c("gene_id", "gene_name"))
sites <- annotate_gene_feature(sites, site_gr, exon_gr, "exon")
intron_hits <- derive_intron_overlaps(
    gene_body_hits, exon_hits,
    unique(as.data.table(gene_table[, c("gene_id", "gene_name")]))
)
sites <- add_overlap_columns(sites, intron_hits, "n_intron_genes", c("gene_id", "gene_name"))
setnames(sites, c("gene_id", "gene_name"), c("intron_gene_id", "intron_gene_name"))
sites[, in_intron := n_intron_genes > 0]
sites <- annotate_gene_feature(sites, site_gr, cds_gr, "cds")
sites <- annotate_gene_feature(sites, site_gr, five_prime_utr_gr, "five_prime_utr")
sites <- annotate_gene_feature(sites, site_gr, three_prime_utr_gr, "three_prime_utr")
sites <- annotate_gene_feature(sites, site_gr, utr_gr, "utr")
sites[, genomic_context := fcase(in_promoter, "promoter", in_gene_body, "gene_body", default = "intergenic")]
sites[, gene_feature_context := fcase(
    in_promoter, "promoter",
    in_exon, "exon",
    in_intron, "intron",
    in_gene_body, "gene_body_unclassified",
    default = "intergenic"
)]

message("Loading and annotating UCSC CpG-island context")
cpg_islands <- read_cpg_islands(opt$CpGIslandAnnotations)
cpg_island_gr <- plyranges::as_granges(data.frame(
    seqnames = as.character(cpg_islands$chrom),
    start = cpg_islands$chromStart + 1L,
    end = cpg_islands$chromEnd,
    strand = "*",
    cpg_island_name = as.character(cpg_islands$name)
))
cpg_island_hits <- summarize_overlaps(site_gr, cpg_island_gr, "cpg_island_name")
sites <- add_overlap_columns(sites, cpg_island_hits, "n_overlapping_cpg_islands", "cpg_island_name")
sites[, in_cpg_island := n_overlapping_cpg_islands > 0]

nearest_cpg_island <- GenomicRanges::nearest(site_gr, cpg_island_gr, ignore.strand = TRUE)
sites[, `:=`(nearest_cpg_island_distance = NA_integer_, nearest_cpg_island_name = NA_character_)]
has_nearest_cpg_island <- !is.na(nearest_cpg_island)
if (any(has_nearest_cpg_island)) {
    queries <- which(has_nearest_cpg_island)
    subjects <- nearest_cpg_island[has_nearest_cpg_island]
    site_positions <- GenomicRanges::start(site_gr)[queries]
    sites[queries, `:=`(
        nearest_cpg_island_distance = as.integer(pmax(
            GenomicRanges::start(cpg_island_gr)[subjects] - site_positions,
            site_positions - GenomicRanges::end(cpg_island_gr)[subjects],
            0L
        )),
        nearest_cpg_island_name = S4Vectors::mcols(cpg_island_gr)$cpg_island_name[subjects]
    )]
}
within_shore_window <- overlaps_any(site_gr, expand_ranges(cpg_island_gr, 2000L))
within_shelf_window <- overlaps_any(site_gr, expand_ranges(cpg_island_gr, 4000L))
sites[, cpg_island_context := fcase(
    in_cpg_island, "island",
    within_shore_window, "shore",
    within_shelf_window, "shelf",
    default = "open_sea"
)]

message("Loading and annotating ENCODE cCRE overlaps")
ccres <- fread(opt$CCREAnnotations, header = FALSE)
if (ncol(ccres) < 6) stop("cCRE annotations must have at least six columns (V1-V6)")
setnames(ccres, names(ccres)[seq_len(6)], c("ccre_chrom", "ccre_begin", "ccre_end", "ccre_v4_id", "ccre_v5_id", "ccre_type"))
ccres[, `:=`(ccre_begin = as.integer(ccre_begin), ccre_end = as.integer(ccre_end))]
if (anyNA(ccres$ccre_begin) || anyNA(ccres$ccre_end) || any(ccres$ccre_begin < 0) || any(ccres$ccre_end <= ccres$ccre_begin)) {
    stop("cCRE annotations must use valid BED start/end coordinates")
}
ccre_gr <- plyranges::as_granges(data.frame(
    seqnames = as.character(ccres$ccre_chrom), start = ccres$ccre_begin + 1L,
    end = ccres$ccre_end, strand = "*", ccre_v4_id = as.character(ccres$ccre_v4_id),
    ccre_v5_id = as.character(ccres$ccre_v5_id), ccre_type = as.character(ccres$ccre_type)
))
ccre_hits <- summarize_overlaps(site_gr, ccre_gr, c("ccre_v4_id", "ccre_v5_id", "ccre_type"))
sites <- add_overlap_columns(sites, ccre_hits, "n_overlapping_ccres", c("ccre_v4_id", "ccre_v5_id", "ccre_type"))
sites[, `:=`(
    in_ccre = n_overlapping_ccres > 0,
    is_enhancer_like = !is.na(ccre_type) & grepl("(^|,)(pELS|dELS)(,|$)", ccre_type),
    is_ctcf_only = !is.na(ccre_type) & grepl("(^|,)CTCF-only(,|$)", ccre_type)
)]

setorder(sites, `#chrom`, begin, end)
sites[, site_index := NULL]
fwrite(sites, output_file, sep = "\t", na = "NA")
message("Wrote annotations for ", nrow(sites), " passing sites: ", output_file)
