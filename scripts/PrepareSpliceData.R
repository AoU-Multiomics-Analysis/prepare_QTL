library(tidyverse)
library(data.table)
library(magrittr)
library(optparse)
library(data.table)
library(rtracklayer)
library(RNOmni)


transform_phenotype <- function(x, rank_normalize){
    if (rank_normalize) {
        return(RankNorm(x))
    }
    as.numeric(scale(x, center = TRUE, scale = TRUE))
}


######## COMMAND LINE ARGUMENTS ############

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
    optparse::make_option(c("--SpliceData"), type="character", default=NULL,
                        help="Bedfile from leafcutter", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--SampleList"), type="character", default=NULL,
                        help="File containing list of samples to run processing on", metavar = "type"),
    optparse::make_option(c("--RankNormalize"), type="character", default="true",
                        help="Deprecated; both INT and scaled BED outputs are always written", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


IntOutputFile <- paste0(opt$OutputPrefix,'.splicing.INT.bed.gz')
ScaledOutputFile <- paste0(opt$OutputPrefix,'.splicing.scaled.bed.gz')
RawOutputFile <- paste0(opt$OutputPrefix,'.splicing.raw.bed.gz')
message(paste0('Writing INT bed file to: ', IntOutputFile))
message(paste0('Writing scaled bed file to: ', ScaledOutputFile))
message(paste0('Writing raw bed file to: ', RawOutputFile))

OutputFileGroups <- paste0(opt$OutputPrefix,'.phenotype_groups.tsv')
message(paste0('Writing to groups file: ',OutputFileGroups))


############### LOAD DATA ###################

message('Loading sample list')
# load sample list data
SampleList <-  readr::read_tsv(opt$SampleList) %>% dplyr::rename('ID' = 1) %>% mutate(ID = as.character(ID))  %>% pull(ID)
message(paste0('Number of samples in SampleList:',SampleList %>% length()))


message('Loading splice data')
SpliceData <-  fread(opt$SpliceData) %>%
    dplyr::select(1,2,3,4,any_of(SampleList))
NumSampleSpliceData <- SpliceData %>% ncol - 4
message(paste0('Number of samples found in SpliceData matching SampleList:', NumSampleSpliceData ))


message('Extracting interval information')
# extracts interval information for splice junctions
SpliceDataTSS <- SpliceData %>% select(1,2,3,4)
SpliceDataValues <- SpliceData %>% dplyr::select(-1,-2,-3,-4)


write_splice_bed <- function(splice_data_values, splice_data_tss, output_file, transform_label, rank_normalize = NULL){
    if (is.null(rank_normalize)) {
        message('Preparing raw splice data')
        SpliceDataNorm <- splice_data_values %>%
            data.frame() %>%
            dplyr::rename_with(~str_remove(.,'X'))
    } else {
        message(paste0('Applying ', transform_label, ' transformation to splice data'))
        # transforms each splice junction across samples
        SpliceDataNorm <- splice_data_values %>%
            t() %>%
            data.frame() %>%
            mutate(across(everything(),~transform_phenotype(., rank_normalize))) %>%
            t() %>%
            data.frame() %>%
            dplyr::rename_with(~str_remove(.,'X'))
    }

    message('Merging normalized splice data and TSS info')
    SpliceDataBed <- bind_cols(splice_data_tss,SpliceDataNorm) %>%
        #arrange(`#chr`,start) %>%
        dplyr::rename('phenotype_id' = 'ID')
    #    mutate(group_id = str_remove(phenotype_id,".*(?=ENSG)"))  %>%
        #arrange(group_id, `#chr`, start, end, phenotype_id)

    message('Extracting phenotype groups')
    #PhenotypeGroups <- SpliceDataBed %>%
            #select(phenotype_id) %>%
            #select(phenotype_id,group_id)

    # Extract ENSG, Sort by CHR, and then ENSG
    SpliceDataBed <- SpliceDataBed %>%
      mutate(ENSG_id = str_extract(phenotype_id, "ENSG[0-9]+"))
    SpliceDataBed <- SpliceDataBed %>%
      arrange(`#chr`, start, end, ENSG_id) %>%
      select(-ENSG_id)

    message(paste0('Writing ', transform_label, ' bedfile to ', output_file))
    SpliceDataBed %>%  fwrite(output_file,sep ='\t')
}

write_splice_bed(SpliceDataValues, SpliceDataTSS, IntOutputFile, 'rank-normalized', TRUE)
write_splice_bed(SpliceDataValues, SpliceDataTSS, ScaledOutputFile, 'scaled', FALSE)
write_splice_bed(SpliceDataValues, SpliceDataTSS, RawOutputFile, 'raw')

#message('Writing phenotype groups')
#SpliceDataBed %>% select(phenotype_id,group_id)%>% fwrite(OutputFileGroups,sep ='\t',col.names = FALSE)
