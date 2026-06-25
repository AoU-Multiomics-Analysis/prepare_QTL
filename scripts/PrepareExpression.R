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

transform_phenotype <- function(x, rank_normalize){
    if (rank_normalize) {
        return(RankNorm(x))
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
 #   optparse::make_option(c("--TPMGCT"), type="character", default=NULL,
                        #help="Parquet or TSV of normalzied protein expression data", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--AnnotationGTF"), type="character", default=NULL,
                        help="GTF file used to TSS locations for each gene", metavar = "type"),
    optparse::make_option(c("--SampleList"), type="character", default=NULL,
                        help="File containing list of samples to run processing on", metavar = "type"),
    optparse::make_option(c("--RankNormalize"), type="character", default="true",
                        help="Deprecated; both INT and scaled BED outputs are always written", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


########### LOAD DATA #####################

message('Loading count data')
CountData <-  fread(opt$CountGCT,skip  =2 ,header = TRUE)

#TPMData <- fread(opt$TPMGCT,skip  =2 ,header = TRUE)

PositionTSS <- extract_TSS_pos(opt$AnnotationGTF)

SampleList <- fread(opt$SampleList,header = FALSE) %>% dplyr::rename('ID' = 1) %>% pull(ID)
nSamples <- SampleList %>% length()
message(paste0('Number of sample in sample list:',nSamples))


IntOutputFile <- paste0(opt$OutputPrefix,'.expression.INT.bed.gz')
ScaledOutputFile <- paste0(opt$OutputPrefix,'.expression.scaled.bed.gz')
RawOutputFile <- paste0(opt$OutputPrefix,'.expression.raw.bed.gz')
message(paste0('Writing INT bed file to: ', IntOutputFile ))
message(paste0('Writing scaled bed file to: ', ScaledOutputFile ))
message(paste0('Writing raw bed file to: ', RawOutputFile ))
############# PROCESS DATA ###########
# transpose read count data such that
# genes are column and rows are samples
message('Transposing data')
CountDataTransposed <- CountData %>%
    dplyr::select(-Description) %>%
    column_to_rownames('Name') %>%
    dplyr::select(any_of(SampleList)) %>%
    t() %>%
    data.frame()

# filter to genes where the count is greater than 6
# in atleast 20% of samples
message('Filtering expression by counts')
CountDataFiltered <- CountDataTransposed %>%
        dplyr::select(where(~ mean(.x > 6) >= 0.2)) %>%
        t() %>%
        data.frame()

message('Performing edgeR TMM normalization')
# Convert filtered count data to DGE list for normalization
DataEdgeR <- edgeR::DGEList(CountDataFiltered)
DataEdgeR <- edgeR::calcNormFactors(DataEdgeR)

message('Computing CPMs')
DataCPM <- edgeR::cpm(DataEdgeR, log=FALSE) %>% data.frame()

write_expression_bed <- function(cpm_data, tss_positions, output_file, transform_label, rank_normalize = NULL, remove_outliers = TRUE){
    message(paste0('Preparing ', transform_label, ' CPM BED'))
    if (is.null(rank_normalize)) {
        NormalizedCPMsMatrix <- cpm_data %>%
                        data.frame() %>%
                        dplyr::rename_with(~str_remove(.,'^X'))
    } else {
        NormalizedCPMsMatrix <- cpm_data %>%
                        t() %>%
                        data.frame() %>%
                        mutate(across(everything(),~transform_phenotype(., rank_normalize))) %>%
                        t() %>%
                        data.frame() %>%
                        dplyr::rename_with(~str_remove(.,'^X'))
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
                inner_join(NormalizedCPMs,by = 'gene_id') %>%
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
        arrange(seqnames,start) %>%
        dplyr::rename('#chr' = 'seqnames') %>%
        fwrite(output_file,sep='\t')
}

write_expression_bed(DataCPM, PositionTSS, IntOutputFile, 'rank-normalized', TRUE)
write_expression_bed(DataCPM, PositionTSS, ScaledOutputFile, 'scaled', FALSE)
write_expression_bed(DataCPM, PositionTSS, RawOutputFile, 'raw', remove_outliers = FALSE)
