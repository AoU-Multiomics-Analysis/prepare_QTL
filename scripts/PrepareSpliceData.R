library(tidyverse)
library(data.table)
library(magrittr)
library(optparse)
library(data.table)
library(rtracklayer)
library(RNOmni)


######## COMMAND LINE ARGUMENTS ############

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
    optparse::make_option(c("--SpliceData"), type="character", default=NULL,
                        help="Bedfile from leafcutter", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--SampleList"), type="character", default=NULL,
                        help="File containing list of samples to run processing on", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


OutputFile <- paste0(opt$OutputPrefix,'.splicing.bed.gz')
message(paste0('Writing to output file: ',OutputFile))

############### LOAD DATA ###################

message('Loading sample list')
# load sample list data
SampleList <-  readr::read_tsv(opt$SampleList) %>% dplyr::rename('ID' = 1) %>% pull(ID)
message(paste0('Number of samples in SampleList:',SampleList %>% length()))


message('Loading splice data')
SpliceData <-  fread(opt$SpliceData) %>% 
    dplyr::select(1,2,3,4,any_of(SampleList))
NumSampleSpliceData <- SpliceData %>% ncol - 4
message(paste0('Number of samples found in SpliceData matching SampleList:', NumSampleSpliceData ))


# extracts interval information for splice junctions
SpliceDataTSS <- SpliceData %>% select(1,2,3,4)


# drops splice interval information and 
# performs RankNorm transformation
SpliceDataNorm <- SpliceData %>% 
    dplyr::select(-1,-2,-3,-4) %>% 
    t() %>% 
    data.frame() %>% 
    mutate(across(everything(),~RankNorm(.))) %>% 
    t() %>% 
    data.frame() %>%
    dplyr::rename_with(~str_remove(.,'X')) 

SpliceDataBed <- bind_cols(SpliceDataTSS,SpliceDataNorm) %>% arrange(seqnames,start)
SpliceDataBed %>% fwrite(OutputFile,sep ='t')

