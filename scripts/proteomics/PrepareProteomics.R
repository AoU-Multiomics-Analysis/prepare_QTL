library(tidyverse)
library(data.table)
library(OlinkAnalyze)
library(biomaRt)
library(optparse)
library(arrow)
library(rtracklayer)
library(RNOmni)
library(WGCNA)


#########  FUNCTIONS ##########

# helper functions that converts ensembl IDs
# to UniProt used in the bedfile annotation.
# This could error out depending on the GTF used
# and how that matches with ensembl version. GENCODE v48
# uses ensembl 114 so thats what im using by default here
check_required_columns <- function(data, required_columns, data_name){
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) {
    stop(paste0(data_name, " is missing required columns: ", paste(missing_columns, collapse = ", ")))
}
}

select_sample_id_column <- function(data, sample_list){
sample_id_columns <- intersect(c("SampleID", "ResearchID"), colnames(data))
if (length(sample_id_columns) == 0) {
    stop("ProteomicData must contain at least one sample ID column: SampleID or ResearchID")
}
sample_list <- as.character(sample_list)
overlap_counts <- sapply(sample_id_columns, function(column) {
    sum(as.character(data[[column]]) %in% sample_list)
})
selected_column <- names(which.max(overlap_counts))
if (overlap_counts[[selected_column]] == 0) {
    stop("No samples in SampleList were found in ProteomicData SampleID or ResearchID columns")
}
message(paste0("Using ", selected_column, " as proteomics sample ID column; matched ", overlap_counts[[selected_column]], " rows"))
selected_column
}

GetUniProtConversion <- function(proteomics_data, ensembl_version = 114){
message(paste0('Using ensembl version ',ensembl_version))
# use biomart to map unprot ids to ensembl gene ids
protein_list <- proteomics_data %>% dplyr::select(UniProt) %>% distinct() %>% pull(UniProt)
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl",version = ensembl_version)

ensembl <- useDataset("hsapiens_gene_ensembl", mart = mart)
conversion_list <- getBM(c("ensembl_gene_id_version","uniprot_gn_id"), "uniprot_gn_id", protein_list, ensembl)
conversion_list <- conversion_list %>% dplyr::rename('gene_id' = 1,'UniProt' = 2)
conversion_list

}


# use rtracklayer to import GTF file and extract TSS locations.
# This should run on the collapsed GTF that has been generated
# but might work with any GTF
extract_TSS_pos <- function(gencode_file){

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
  #TODO look around if there is a package recognizing delimiter in dataset
    optparse::make_option(c("--ProteomicData"), type="character", default=NULL,
                        help="Parquet or TSV of normalzied protein expression data", metavar = "type"),
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

filepath <- opt$ProteomicData
############### LOAD DATA ###################
# load in proteomics data based on type of file
if (grepl("\\.tsv(\\.gz)?$", filepath)) {
    message("Reading as TSV using fread")
    ProteomicsData <-  readr::read_tsv(filepath)
  } else if (grepl("\\.parquet$", filepath)) {
    message("Reading as Parquet using arrow")
    ProteomicsData <- arrow::read_parquet(filepath)
  } else {
    stop("Unsupported file extension. Use .tsv, .tsv.gz, or .parquet")
  }

check_required_columns(ProteomicsData, c("PCNormalizedNPX", "UniProt"), "ProteomicData")

message('Loading sample list')
# load sample list data
SampleList <-  readr::read_tsv(opt$SampleList) %>% dplyr::rename('ID' = 1) %>% mutate(ID = as.character(ID)) %>% pull(ID)
SampleIDColumn <- select_sample_id_column(ProteomicsData, SampleList)

message('Creating UniProt ensembl id conversion')
# create table containing UniProt and ensembl id conversions
UniProtConversion <- GetUniProtConversion(ProteomicsData)
check_required_columns(UniProtConversion, c("gene_id", "UniProt"), "UniProt conversion")

message('Extracting TSS locations from GTF')
# extract TSS locations from input GTF
TSS_position <- extract_TSS_pos(opt$AnnotationGTF)
check_required_columns(TSS_position, "gene_id", "TSS position table")

message('Merging UniProt conversion table and TSS locations')
UniProtTSSLocations <- UniProtConversion %>% left_join(TSS_position,by = 'gene_id')

IntOutputFile <- paste0(opt$OutputPrefix,'.protein.INT.bed.gz')
ScaledOutputFile <- paste0(opt$OutputPrefix,'.protein.scaled.bed.gz')
RawOutputFile <- paste0(opt$OutputPrefix,'.protein.raw.bed.gz')
message(paste0('Writing INT bed file to ', IntOutputFile))
message(paste0('Writing scaled bed file to ', ScaledOutputFile))
message(paste0('Writing raw bed file to ', RawOutputFile))
############ BEGIN SCRIPT ##########


message('Converting data to wide format')
# converts data to wide format such that
# each column is a protein and each row is the
# quantification in an individual.
ProteomicsDataWideRaw <- ProteomicsData %>%
    filter(as.character(.data[[SampleIDColumn]]) %in% SampleList) %>%
    #group_by(ResearchID,OlinkID) %>%
    #filter(row_number() == 1) %>%
    #ungroup() %>%
    transmute(SampleIDForQTL = as.character(.data[[SampleIDColumn]]), PCNormalizedNPX, UniProt) %>%
    pivot_wider(names_from = UniProt,values_from =PCNormalizedNPX )  %>%
    column_to_rownames('SampleIDForQTL')


write_proteomics_bed <- function(proteomics_data_wide, uniprot_tss_locations, output_file, transform_label, rank_normalize = NULL, remove_outliers = TRUE){
    if (is.null(rank_normalize)) {
        message('Preparing raw proteomics data')
        ProteomicsDataWide <- proteomics_data_wide
    } else {
        message(paste0('Applying ', transform_label, ' transformation to proteomics data'))
        ProteomicsDataWide <- proteomics_data_wide %>%
            mutate(across(everything(),~transform_phenotype(., rank_normalize)))
    }

    message('Merging data with TSS locations')
    # converts normalized proteomic data to
    # bed format where each row is a gene/protein
    # and columns contain interval TSS location, molecular trait id
    # and other columns are participant quantifications
    ProteomicsFeatureMatrix <- ProteomicsDataWide %>%
        t() %>%
        data.frame() %>%
        dplyr::rename_with(~str_remove(.,'^X'))

    if (remove_outliers) {
        ProteomicsFeatureMatrix <- remove_connectivity_outliers(ProteomicsFeatureMatrix, output_file, transform_label)
    } else {
        message(paste0('Skipping connectivity outlier removal for ', transform_label, ' data; keeping ', ncol(ProteomicsFeatureMatrix), ' samples'))
    }

    ProteomicsBed <- ProteomicsFeatureMatrix %>%
        rownames_to_column('UniProt') %>%
        left_join(uniprot_tss_locations,by = 'UniProt') %>%
        dplyr::select(seqnames,start, end, UniProt, gene_id, everything()) %>%
        filter(!is.na(seqnames)) %>%
        dplyr::rename_with(~str_remove(.,'^X')) %>%
        mutate(gene_id = paste0(UniProt, '_', gene_id)) %>%
        dplyr::select(-UniProt) %>%
        arrange(seqnames,start) %>%
        dplyr::rename('#chr'='seqnames')

    nproteins <- ProteomicsBed %>% nrow
    message(paste0(nproteins, ' proteins found'))
    message(paste0('Writing ', transform_label, ' data to ', output_file))
    ProteomicsBed %>% fwrite(output_file,sep ='\t')
}


write_proteomics_bed(ProteomicsDataWideRaw, UniProtTSSLocations, IntOutputFile, 'rank-normalized', TRUE)
write_proteomics_bed(ProteomicsDataWideRaw, UniProtTSSLocations, ScaledOutputFile, 'scaled', FALSE)
write_proteomics_bed(ProteomicsDataWideRaw, UniProtTSSLocations, RawOutputFile, 'raw', remove_outliers = FALSE)
