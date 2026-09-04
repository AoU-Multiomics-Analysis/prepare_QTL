library(tidyverse)
library(data.table)
library(magrittr)
library(biomaRt)
library(optparse)
library(data.table)
library(rtracklayer)
library(RNOmni)
library(edgeR)
library(WGCNA)


read_expression_sample_list <- function(path, expression_samples, expression_label) {
    sample_table <- readr::read_tsv(path, col_names = FALSE,
        col_types = readr::cols(.default = readr::col_character()),
        na = character(), trim_ws = TRUE, skip_empty_rows = FALSE,
        show_col_types = FALSE, progress = FALSE)
    if (ncol(sample_table) == 0L || nrow(sample_table) == 0L) {
        stop('Sample list is empty', call. = FALSE)
    }
    if (ncol(sample_table) != 1L || nrow(readr::problems(sample_table)) > 0L) {
        stop('Sample list must contain exactly one column', call. = FALSE)
    }
    samples <- dplyr::pull(sample_table, 1)
    if (anyNA(samples) || any(!nzchar(samples))) {
        stop('Sample list contains a blank sample ID', call. = FALSE)
    }
    first <- samples[[1L]]
    if (tolower(first) %in% c('research_id', 'sample_id', 'id') &&
        !(first %in% expression_samples)) {
        message(paste0("Sample list: removed header '", first, "'"))
        samples <- samples[-1L]
    }
    if (length(samples) == 0L) {
        stop('Sample list contains no sample IDs after removing the header', call. = FALSE)
    }
    if (anyDuplicated(samples)) {
        stop('Sample list contains duplicate sample IDs', call. = FALSE)
    }
    missing_samples <- setdiff(samples, expression_samples)
    if (length(missing_samples) > 0L) {
        stop(paste0('Samples missing from ', expression_label, ': ',
            paste(missing_samples, collapse = ', ')), call. = FALSE)
    }
    message(paste0('Number of sample in sample list:', length(samples)))
    samples
}


# use rtracklayer to import GTF file and extract TSS locations.
# This should run on the collapsed GTF that has been generated
# but might work with any GTF
extract_TSS_pos <- function(gencode_file){
message('Extracting TSS locations')
message('Loading GTF')
gencode_GTF <- rtracklayer::import(gencode_file) %>% data.frame()

message('Extracting GTF')
# map TSS locations based on strand
TSS_locations <- gencode_GTF  %>%
    filter(type == 'gene'  ) %>%
    mutate(TSS = case_when(strand == '+' ~ start,TRUE ~ end)) %>%
    dplyr::select(gene_id,TSS,seqnames) %>%
    mutate(start = TSS -1,end = TSS) %>%
    dplyr::select(gene_id,start,end,seqnames)
TSS_locations
}

transform_phenotype <- function(x, rank_normalize, log2_transform = FALSE){
    if (rank_normalize) {
        return(RankNorm(x))
    }
    if (log2_transform) {
        x <- log2(x + 1)
    }
    as.numeric(scale(x, center = TRUE, scale = TRUE))
}

remove_connectivity_outliers <- function(phenotype_matrix, output_file, transform_label){
    outliers_file <- str_replace(output_file, "\\.bed\\.gz$", ".connectivity_outliers.tsv")
    message(paste0('Computing connectivity outliers for ', transform_label, ' data'))

    n_samples_before <- ncol(phenotype_matrix)
    empty_outliers <- tibble(SampleID = character(), Z_score = numeric())
    if (ncol(phenotype_matrix) < 3 || nrow(phenotype_matrix) < 2) {
        message('Not enough data to compute connectivity outliers; keeping all samples')
        message(paste0('Connectivity outlier removal for ', transform_label, ' data: removed 0 of ', n_samples_before, ' samples; ', n_samples_before, ' samples remain'))
        empty_outliers %>% write_tsv(outliers_file)
        return(phenotype_matrix)
    }

    phenotype_matrix <- as.data.frame(phenotype_matrix, check.names = FALSE)
    norm_adj <- 0.5 + 0.5 * WGCNA::bicor(phenotype_matrix, use = "pairwise.complete.obs")
    norm_adj[is.na(norm_adj)] <- 0

    net_summary <- WGCNA::fundamentalNetworkConcepts(norm_adj)
    net_connectivity <- net_summary$Connectivity
    connectivity_sd <- sd(net_connectivity, na.rm = TRUE)

    if (is.na(connectivity_sd) || connectivity_sd == 0) {
        message('Connectivity scores have zero or undefined variance; keeping all samples')
        message(paste0('Connectivity outlier removal for ', transform_label, ' data: removed 0 of ', n_samples_before, ' samples; ', n_samples_before, ' samples remain'))
        empty_outliers %>% write_tsv(outliers_file)
        return(phenotype_matrix)
    }

    connectivity_zscore <- ((net_connectivity - mean(net_connectivity, na.rm = TRUE)) / connectivity_sd) %>%
        data.frame() %>%
        dplyr::rename('Z_score' = 1) %>%
        rownames_to_column('SampleID')

    connectivity_zscore_outliers <- connectivity_zscore %>% filter(Z_score < -3)
    n_samples_removed <- nrow(connectivity_zscore_outliers)
    n_samples_after <- n_samples_before - n_samples_removed
    message(paste0('Connectivity outlier removal for ', transform_label, ' data: removed ', n_samples_removed, ' of ', n_samples_before, ' samples; ', n_samples_after, ' samples remain'))
    connectivity_zscore_outliers %>% write_tsv(outliers_file)

    kept_samples <- setdiff(colnames(phenotype_matrix), connectivity_zscore_outliers$SampleID)
    phenotype_matrix[, kept_samples, drop = FALSE]
}



######## COMMAND LINE ARGUMENTS ############

option_list <- list(
    optparse::make_option(c("--CountGCT"), type="character", default=NULL,
                        help="Parquet or TSV of normalzied protein expression data", metavar = "type"),
    optparse::make_option(c("--CpmBed"), type="character", default=NULL,
                        help="Coordinate-preserving BED of pre-normalized linear CPM values", metavar = "type"),
    optparse::make_option(c("--Log2CpmBed"), type="character", default=NULL,
                        help="Coordinate-preserving BED of pre-normalized log2 CPM values", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--AnnotationGTF"), type="character", default=NULL,
                        help="GTF file used to TSS locations for each gene", metavar = "type"),
    optparse::make_option(c("--SampleList"), type="character", default=NULL,
                        help="One sample ID per line; optional research_id, sample_id, or ID header", metavar = "type"),
    optparse::make_option(c("--RankNormalize"), type="character", default="true",
                        help="Deprecated; both INT and scaled BED outputs are always written", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


has_count_gct <- !is.null(opt$CountGCT) && nzchar(opt$CountGCT)
has_cpm_bed <- !is.null(opt$CpmBed) && nzchar(opt$CpmBed)
has_log2_cpm_bed <- !is.null(opt$Log2CpmBed) && nzchar(opt$Log2CpmBed)

if (sum(c(has_count_gct, has_cpm_bed, has_log2_cpm_bed)) != 1L) {
    stop('Provide exactly one of --CountGCT, --CpmBed, or --Log2CpmBed')
}
if (has_count_gct && (is.null(opt$AnnotationGTF) || !nzchar(opt$AnnotationGTF))) {
    stop('--AnnotationGTF is required with --CountGCT')
}
if (!has_count_gct && !is.null(opt$AnnotationGTF) && nzchar(opt$AnnotationGTF)) {
    bed_option <- if (has_cpm_bed) '--CpmBed' else '--Log2CpmBed'
    stop(paste('--AnnotationGTF cannot be used with', bed_option))
}


########### LOAD DATA #####################

IntOutputFile <- paste0(opt$OutputPrefix,'.expression.INT.bed.gz')
ScaledOutputFile <- paste0(opt$OutputPrefix,'.expression.scaled.bed.gz')
RawOutputFile <- paste0(opt$OutputPrefix,'.expression.raw.bed.gz')
message(paste0('Writing INT bed file to: ', IntOutputFile ))
message(paste0('Writing scaled bed file to: ', ScaledOutputFile ))
message(paste0('Writing raw bed file to: ', RawOutputFile ))
############# PROCESS DATA ###########
if (has_count_gct) {
    message('Loading count data')
    CountData <- fread(opt$CountGCT, skip = 2, header = TRUE)
    SampleList <- read_expression_sample_list(opt$SampleList,
        setdiff(names(CountData), c('Name', 'Description')), 'count GCT')
    PositionTSS <- extract_TSS_pos(opt$AnnotationGTF)

    # transpose read count data such that
    # genes are column and rows are samples
    message('Transposing data')
    CountDataTransposed <- CountData %>%
        dplyr::select(-Description) %>%
        column_to_rownames('Name') %>%
        dplyr::select(any_of(SampleList)) %>%
        t() %>%
        data.frame(check.names = FALSE)

    # filter to genes where the count is greater than 6
    # in atleast 20% of samples
    message('Filtering expression by counts')
    CountDataFiltered <- CountDataTransposed %>%
            dplyr::select(where(~ mean(.x > 6) >= 0.2)) %>%
            t() %>%
            data.frame(check.names = FALSE)

    message('Performing edgeR TMM normalization')
    # Convert filtered count data to DGE list for normalization
    DataEdgeR <- edgeR::DGEList(CountDataFiltered)
    DataEdgeR <- edgeR::calcNormFactors(DataEdgeR)

    message('Computing CPMs')
    DataCPM <- edgeR::cpm(DataEdgeR, log = FALSE) %>% data.frame(check.names = FALSE)
} else {
    bed_label <- if (has_cpm_bed) 'CPM BED' else 'Log2 CPM BED'
    message(paste('Loading pre-normalized', bed_label))
    ExpressionBed <- fread(
        if (has_cpm_bed) opt$CpmBed else opt$Log2CpmBed,
        header = TRUE,
        check.names = FALSE,
        colClasses = list(character = c('#chr', 'gene_id'))
    )
    metadata_columns <- c('#chr', 'start', 'end', 'gene_id')
    if (!identical(names(ExpressionBed)[seq_along(metadata_columns)], metadata_columns)) {
        stop(paste(bed_label, 'must start with #chr, start, end, and gene_id columns'))
    }
    if (anyDuplicated(ExpressionBed$gene_id)) {
        stop(paste(bed_label, 'gene_id values must be unique'))
    }
    SampleList <- read_expression_sample_list(opt$SampleList,
        names(ExpressionBed)[-(seq_along(metadata_columns))], bed_label)

    DataCPM <- ExpressionBed %>%
        dplyr::select(gene_id, all_of(SampleList)) %>%
        column_to_rownames('gene_id') %>%
        data.frame(check.names = FALSE)
    if (!all(vapply(DataCPM, is.numeric, logical(1))) || any(!is.finite(as.matrix(DataCPM)))) {
        stop(paste(bed_label, 'sample values must be finite numeric values'))
    }
    if (has_cpm_bed && any(as.matrix(DataCPM) < 0)) {
        negative_genes <- rownames(DataCPM)[rowSums(as.matrix(DataCPM) < 0) > 0]
        stop(paste0(
            'CPM BED contains negative values; linear CPM must be nonnegative. ',
            'Values were not clipped. Example gene IDs: ',
            paste(head(negative_genes, 10), collapse = ', ')
        ))
    }
    zero_variance_genes <- DataCPM %>%
        rownames_to_column('gene_id') %>%
        rowwise() %>%
        mutate(.gene_sd = sd(c_across(-gene_id))) %>%
        ungroup() %>%
        filter(is.na(.gene_sd) | .gene_sd == 0) %>%
        pull(gene_id)
    if (length(zero_variance_genes) > 0) {
        displayed_genes <- head(zero_variance_genes, 10)
        stop(paste0(
            bed_label, ' contains genes with zero variance: ',
            paste(displayed_genes, collapse = ', ')
        ))
    }

    PositionTSS <- ExpressionBed %>%
        dplyr::select(all_of(metadata_columns)) %>%
        dplyr::rename(seqnames = '#chr') %>%
        dplyr::mutate(.input_order = dplyr::row_number())
    message('Skipping count filtering, TMM normalization, CPM calculation, and GTF mapping')
}

write_expression_bed <- function(cpm_data, tss_positions, output_file, transform_label, rank_normalize = NULL, log2_transform = FALSE, remove_outliers = TRUE){
    message(paste0('Preparing ', transform_label, ' CPM BED'))
    if (is.null(rank_normalize)) {
        NormalizedCPMsMatrix <- cpm_data %>%
                        data.frame(check.names = FALSE)
    } else {
        if (log2_transform) {
            message('Applying log2(CPM + 1) before centering and scaling')
        }
        NormalizedCPMsMatrix <- cpm_data %>%
                        t() %>%
                        data.frame(check.names = FALSE) %>%
                        mutate(across(everything(),~transform_phenotype(., rank_normalize, log2_transform))) %>%
                        t() %>%
                        data.frame(check.names = FALSE)
    }

    if (remove_outliers) {
        NormalizedCPMsMatrix <- remove_connectivity_outliers(NormalizedCPMsMatrix, output_file, transform_label)
    } else {
        message(paste0('Skipping connectivity outlier removal for ', transform_label, ' data; keeping ', ncol(NormalizedCPMsMatrix), ' samples'))
    }

    NormalizedCPMs <- NormalizedCPMsMatrix %>% rownames_to_column('gene_id')

    LengthNormalziedCPMS <- NormalizedCPMs %>% nrow
    message(paste0('Number of genes found: ',LengthNormalziedCPMS))

    message('Merging quantifications with TSS locations')
    BedNormalizedCPMs <- tss_positions %>%
                inner_join(NormalizedCPMs,by = 'gene_id')
    if ('.input_order' %in% names(BedNormalizedCPMs)) {
        BedNormalizedCPMs <- BedNormalizedCPMs %>%
            arrange(.input_order) %>%
            dplyr::select(-.input_order)
    } else {
        BedNormalizedCPMs <- BedNormalizedCPMs %>% arrange(seqnames, start)
    }
    BedNormalizedCPMs <- BedNormalizedCPMs %>%
        dplyr::select(seqnames,start,end,gene_id,everything())

    LengthBedCPMs <- BedNormalizedCPMs %>% nrow
    message(paste0('Number of genes found after merge: ',LengthBedCPMs))

    LostGenes <- LengthNormalziedCPMS  - LengthBedCPMs
    message(paste0('Gene lost after merging: ',LostGenes))

    if (LostGenes > 0){
    message('Please check GENCODE version since genes are lost in merging process')

    }

    message(paste0('Writing ', transform_label, ' data to ', output_file))
    BedNormalizedCPMs %>%
        dplyr::rename('#chr' = 'seqnames') %>%
        fwrite(output_file,sep='\t')
}

write_expression_bed(DataCPM, PositionTSS, IntOutputFile, 'rank-normalized', TRUE)
write_expression_bed(DataCPM, PositionTSS, ScaledOutputFile, 'scaled', FALSE, log2_transform = !has_log2_cpm_bed)
write_expression_bed(DataCPM, PositionTSS, RawOutputFile, 'raw', remove_outliers = FALSE)
