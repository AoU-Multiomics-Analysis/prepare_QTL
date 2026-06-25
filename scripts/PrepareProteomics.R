library(tidyverse)
library(data.table)
library(magrittr)
library(biomaRt)
library(optparse)
library(arrow)
library(rtracklayer)
library(RNOmni)


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

GetUniProtConversion <- function(proteomics_data, ensembl_version = 114){
message(paste0('Using ensembl version ',ensembl_version))
# use biomart to map unprot ids to ensembl gene ids
protein_list <- proteomics_data %>% dplyr::select(UniProt) %>% distinct() %>% pull(UniProt)
mart <- useEnsembl("ensembl","hsapiens_gene_ensembl",version = ensembl_version)

ensembl <- useDataset("hsapiens_gene_ensembl", mart = mart)
conversion_list <- getBM(c("ensembl_gene_id_version","uniprot_gn_id"), "uniprot_gn_id", protein_list, ensembl) 
conversion_list %<>% dplyr::rename('gene_id' = 1,'UniProt' = 2)
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

parse_bool <- function(x){
    x <- tolower(as.character(x))
    if (x %in% c("true", "t", "1", "yes", "y")) {
        return(TRUE)
    }
    if (x %in% c("false", "f", "0", "no", "n")) {
        return(FALSE)
    }
    stop("RankNormalize must be true or false")
}

transform_phenotype <- function(x, rank_normalize){
    if (rank_normalize) {
        return(RankNorm(x))
    }
    as.numeric(scale(x, center = TRUE, scale = TRUE))
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
                        help="Apply rank-based inverse normal transformation before writing BED output", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
RankNormalize <- parse_bool(opt$RankNormalize)

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

check_required_columns(ProteomicsData, c("SampleID", "ResearchID", "PCNormalizedNPX", "UniProt"), "ProteomicData")

message('Loading sample list')
# load sample list data
SampleList <-  readr::read_tsv(opt$SampleList) %>% dplyr::rename('ID' = 1) %>% pull(ID)

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

OutputFile <- paste0(opt$OutputPrefix,'.protein.bed.gz')
message(paste0('Writing bed file to ',OutputFile))
############ BEGIN SCRIPT ##########


if (RankNormalize) {
    message('Converting data to wide format and applying rank-based inverse normal transformation')
} else {
    message('Converting data to wide format and centering/scaling')
}
# converts data to wide format such that 
# each column is a protein and each row is the 
# quantification in an individual and then transforms
# each protein across samples.
ProteomicsDataWide <- ProteomicsData %>% 
    filter(SampleID %in% SampleList) %>% 
    #group_by(ResearchID,OlinkID) %>% 
    #filter(row_number() == 1) %>% 
    #ungroup() %>% 
    dplyr::select(ResearchID,PCNormalizedNPX,UniProt) %>% 
    pivot_wider(names_from = UniProt,values_from =PCNormalizedNPX )  %>% 
    column_to_rownames('ResearchID') %>% 
    mutate(across(everything(),~transform_phenotype(., RankNormalize)))


message('Merging data with TSS locations')
# converts normalized proteomic data to 
# bed format where each row is a gene/protein 
# and columns contain interval TSS location, molecular trait id 
# and other columns are  participant  quantifications 
ProteomicsBed <- ProteomicsDataWide %>% 
    t() %>% 
    data.frame() %>%
    dplyr::rename_with(~str_remove(.,'X')) %>% 
    rownames_to_column('UniProt') %>%
    left_join(UniProtTSSLocations,by = 'UniProt') %>% 
    dplyr::select(seqnames,start, end, UniProt, gene_id, everything()) %>% 
    filter(!is.na(seqnames)) %>% 
    dplyr::rename_with(~str_remove(.,'X')) %>% 
    mutate(gene_id = paste0(UniProt, '_', gene_id)) %>% 
    dplyr::select(-UniProt) %>% 
    arrange(seqnames,start) %>% 
    dplyr::rename('#chr'='seqnames')

nproteins <- ProteomicsBed %>% nrow
message('Writing data to output')
message(paste0(nproteins, ' found'))

# write data to output
ProteomicsBed %>% fwrite(OutputFile,sep ='\t')
