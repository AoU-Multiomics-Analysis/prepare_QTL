library(tidyverse)
library(data.table)
library(OlinkAnalyze)
library(magrittr)
library(patchwork)
library(readxl)

############### FUNCTIONS ################
# wrapper function to run median normalization using reference medians
wrap_olink_norm <- function(olink_df,reference_medians){
normed <- olink_normalization(df1 = olink_df %>% as_tibble(),
                    overlapping_samples_df1 = unique(olink_df$SampleID),
                    reference_medians = reference_medians)
normed    
}

########### LOAD DATA #######################

olink_manifest <- 'gs://prod-drc-broad/aou_proteomics/final_10K_20250709/manifest_aou_proteomics_20250709.tsv'
system(paste0('gsutil cp ',olink_manifest,' .'))


olink_10k_set <- '/home/jupyter/multiomics_qtls/edit/olink_10k/parquet'


# load in 10k AoU set
olink_10k_df <-  list.files(path = olink_10k_set,
                    pattern = "parquet$",
                    full.names = TRUE) %>% 
         lapply(OlinkAnalyze::read_NPX)  %>% 
         dplyr::bind_rows()

olink_10k_df_assay_only <- olink_10k_df %>% 
                filter(AssayType == 'assay'& !str_detect(SampleID, "CONTROL_SAMPLE|Control")) 


######### NORMALIZE DATA ################

# use just the first plate by date for reference medians 
reference_medians <- olink_10k_df_assay_only %>% 
                        filter(PlateID == '000171002612_A1_01-17-2024_12-21-52')  %>% 
                        group_by(OlinkID) %>% 
                        summarize(Reference_NPX = median(NPX,na.rm = T))

# adjust all plates by reference median 
normed_data <- olink_10k_df_assay_only %>% 
                        mutate(temp_plate = PlateID) %>% 
                        group_by(temp_plate) %>% 
                        group_modify(~wrap_olink_norm(.,reference_medians)) %>% 
                        select(-temp_plate) %>% 
                        ungroup()

normed_data %>% write_tsv('Olink_10k_df_median_normalized.tsv')


###### FILTER DATA #################
LongOlinkValues <- OlinkMedianNormalized %>% 
    filter(!str_detect(SampleID, 'CONT|Control')  & AssayType == 'assay') %>% 
    mutate(SampleID = str_remove(SampleID,'_.*')) %>% 
    select(SampleID,PlateID,NPX,UniProt)   %>% 
    mutate(joint_id = paste0(SampleID,'_',PlateID))  

valid_sample_plates <- bind_rows(
    LongOlinkValues %>% 
        distinct(SampleID,PlateID) %>% 
        group_by(SampleID) %>% 
        filter(dplyr::n() > 1) %>% 
        filter(row_number() == 1),
    LongOlinkValues %>% 
        distinct(SampleID,PlateID) %>% 
        group_by(SampleID) %>% 
        filter(dplyr::n()  == 1) 
    ) %>% 
    mutate(joint_id = paste0(SampleID,'_',PlateID)) %>% 
    distinct(joint_id)

missing_values <- LongOlinkValues %>% 
    filter(is.na(NPX)) %>% 
    distinct(SampleID) %>% 
    pull(SampleID)

FilteredLongValues <- LongOlinkValues %>% 
    filter(!Sample %in%  missing_values$SampleID) %>% 
    filter(joint_id %in% valid_sample_plates$joint_id) %>% 
    filter(UniProt != 'P32455' & UniProt != 'Q02750') %>% 
    dplyr::rename('PCNormalizedNPX' = 'NPX','ResearchID' = 1) %>% 
    select(ResearchID,UniProt,PCNormalizedNPX)

FilteredLongValues %>% write_tsv('olink_npx_values.tsv.gz')

