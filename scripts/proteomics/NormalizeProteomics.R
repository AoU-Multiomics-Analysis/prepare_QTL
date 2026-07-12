library(tidyverse)
library(OlinkAnalyze)
library(optparse)


############### FUNCTIONS ################

check_required_columns <- function(data, required_columns, data_name){
missing_columns <- setdiff(required_columns, colnames(data))
if (length(missing_columns) > 0) {
    stop(paste0(data_name, " is missing required columns: ", paste(missing_columns, collapse = ", ")))
}
}

# wrapper function to run median normalization using reference medians
wrap_olink_norm <- function(olink_df, reference_medians){
normed <- olink_normalization(df1 = olink_df %>% as_tibble(),
                    overlapping_samples_df1 = unique(olink_df$SampleID),
                    reference_medians = reference_medians)
normed
}


######## COMMAND LINE ARGUMENTS ############

option_list <- list(
    optparse::make_option(c("--OlinkDataDir"), type="character",
                        default="/home/jupyter/multiomics_qtls/edit/olink_10k/parquet",
                        help="Directory containing Olink NPX parquet files", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default="olink",
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--OutputDir"), type="character", default=".",
                        help="Directory for output files", metavar = "type"),
    optparse::make_option(c("--ReferencePlate"), type="character",
                        default="000171002612_A1_01-17-2024_12-21-52",
                        help="PlateID used to calculate reference medians", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

dir.create(opt$OutputDir, recursive = TRUE, showWarnings = FALSE)


########### LOAD DATA #######################

message("Loading Olink NPX parquet files")
olink_files <- list.files(path = opt$OlinkDataDir,
                    pattern = "parquet$",
                    full.names = TRUE)

if (length(olink_files) == 0) {
    stop(paste0("No parquet files found in OlinkDataDir: ", opt$OlinkDataDir))
}

# load in Olink data
olink_df <- olink_files %>%
         lapply(OlinkAnalyze::read_NPX) %>%
         dplyr::bind_rows()

check_required_columns(olink_df, c("SampleID", "PlateID", "AssayType", "OlinkID", "NPX", "UniProt"), "Olink data")

olink_df_assay_only <- olink_df %>%
                filter(AssayType == "assay" & !str_detect(SampleID, "CONTROL_SAMPLE|Control|CONT"))


######### NORMALIZE DATA ################

message("Calculating reference medians")
reference_plate_data <- olink_df_assay_only %>%
                        filter(PlateID == opt$ReferencePlate)

if (nrow(reference_plate_data) == 0) {
    stop(paste0("ReferencePlate was not found in the Olink data: ", opt$ReferencePlate))
}

reference_medians <- reference_plate_data %>%
                        group_by(OlinkID) %>%
                        summarize(Reference_NPX = median(NPX, na.rm = TRUE), .groups = "drop")

message("Normalizing Olink data by plate")
normed_data <- olink_df_assay_only %>%
                        mutate(temp_plate = PlateID) %>%
                        group_by(temp_plate) %>%
                        group_modify(~wrap_olink_norm(., reference_medians)) %>%
                        select(-temp_plate) %>%
                        ungroup()

median_normalized_output <- file.path(opt$OutputDir, paste0(opt$OutputPrefix, "_median_normalized.tsv.gz"))
message(paste0("Writing median-normalized data to ", median_normalized_output))
normed_data %>% write_tsv(median_normalized_output)


###### FILTER DATA #################
warn_fail_assays <- normed_data %>% 
    filter(!str_detect(SampleID, 'CONT|Control')  & AssayType == 'assay') %>% 
    filter(SampleQC != 'PASS' | AssayQC != 'PASS')  

assay_fail_list <- warn_fail_assays %>% 
    filter(AssayQC != 'PASS') %>% 
    distinct(Assay)
sample_remove_list <- warn_fail_assays %>% 
    filter(!Assay %in% assay_fail_list$Assay) %>% 
    filter(SampleQC != 'PASS') %>% 
    distinct(SampleID) 
n_assays_removed <- assay_fail_list  %>% nrow
n_samples_removed <- sample_remove_list  %>% nrow


message('Removing',n_assays_removed,'assays')
message('Removing',n_samples_removed,'samples')



message("Filtering normalized Olink data")
LongOlinkValues <- normed_data %>%
    filter(!str_detect(SampleID, "CONTROL_SAMPLE|Control|CONT") & AssayType == "assay") %>%
    filter(!Assay %in% assay_fail_list$Assay & !SampleID %in% sample_remove_list$SampleID) %>% 
    mutate(SampleID = str_remove(SampleID, "_.*")) %>%
    select(SampleID, PlateID, NPX, UniProt) %>%
    mutate(joint_id = paste0(SampleID, "_", PlateID))

valid_sample_plates <- bind_rows(
    LongOlinkValues %>%
        distinct(SampleID, PlateID) %>%
        group_by(SampleID) %>%
        filter(dplyr::n() > 1) %>%
        filter(row_number() == 1),
    LongOlinkValues %>%
        distinct(SampleID, PlateID) %>%
        group_by(SampleID) %>%
        filter(dplyr::n() == 1)
    ) %>%
    mutate(joint_id = paste0(SampleID, "_", PlateID)) %>%
    distinct(joint_id)

missing_samples <- LongOlinkValues %>%
    filter(is.na(NPX)) %>%
    distinct(SampleID) %>%
    pull(SampleID)

FilteredLongValues <- LongOlinkValues %>%
    filter(!SampleID %in% missing_samples) %>%
    filter(joint_id %in% valid_sample_plates$joint_id) %>%
    filter(!UniProt %in% c("P32455", "Q02750")) %>%
    dplyr::rename("PCNormalizedNPX" = "NPX", "ResearchID" = "SampleID") %>%
    select(ResearchID, UniProt, PCNormalizedNPX)

filtered_output <- file.path(opt$OutputDir, paste0(opt$OutputPrefix, "_npx_values.tsv.gz"))
message(paste0("Writing filtered proteomics data to ", filtered_output))
FilteredLongValues %>% write_tsv(filtered_output)
