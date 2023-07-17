# Script to preprocess metabolomics data
# Generates a file with preprocesed data and a file with metabolite analyzed

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath")))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# input files
metabolomics_input <- "input/mayo_metabolomics.xlsx"
metabolomics_input_checksum <- '_check_sum_of_input_'
metabolomics_metadata <- "input/mayo_metadata.xlsx"

# output files
met_pdata <- 'results/tmp_mayo_metabolomics_processed_data.xlsx'
metabolites_analyzed <- 'results/supplementary_table_X_metabolites_analyzed.xlsx'

# libraries
library(maplet) # omics analysis 

library(tidyverse) # tidy
library(magrittr) # %<>%

library(openxlsx) # excel

source("custom_functions.R") # customized functions

# result directory - creates if does not exist
dir.create("results/", showWarnings = F, recursive = T)

# loading and preprocessing metabolomics data using our default protocol
D  <- 
  # load data and annotations
  mt_load_metabolon_v1(file=metabolomics_input, sheet= "OrigScale") %>%
  # make sure file is unchanged
  mt_load_checksum(file = metabolomics_input, checksum = metabolomics_input_checksum) %>%
  # add annotations
  mt_anno_xls(file = metabolomics_metadata, sheet = 1, anno_type = "samples", data_id_col = "CLIENT IDENTIFIER", anno_id_col = "SampleID_TUBE") %>%
  # create standard variable name SubjectI
  mt_anno_mutate(anno_type ='samples', col_name = 'SubjectID', term=CLIENT.IDENTIFIER) %>% # copy from projid
  # set zeros to NA
  mt_pre_zero_to_na() %>%
  # filter samples with 100% missing values
  mt_pre_filter_missingness(samp_max=0.999999999999) %>% # 99.9..% because the function uses a <= operator
  # filter metabolites with >25% missing values
  mt_pre_filter_missingness(feat_max=0.25) %>%
  # quotient normalization
  mt_pre_norm_quot() %>%
  # log2 transform data
  mt_pre_trans_log() %>%
  # KNN imputation
  mt_pre_impute_knn() %>%
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>% 
  # kNN imputation
  mt_pre_impute_knn() %>%
  # modify the PA name
  mt_anno_mutate(anno_type = 'samples', col_name = 'Diagnosis', 
                 term=case_when(Diagnosis=='Pathologic Aging' ~ 'PA', TRUE ~ Diagnosis)) %>%
  {.}
# naming the assay with sample ids
colnames(D) <- D$SubjectID
# write out the processed data
mt_write_se_xls(D, file=met_pdata)

# write out metabolites i.e supplementary table 2
# metabolites analyzed
mets <- D %>% rowData() %>% data.frame() %>%
  select(-PATHWAY_SORTORDER, - BIOCHEMICAL) %>% dplyr::rename(HMDB=HMDb) %>%
  select(name, everything())
out <- 'metabolites_analyzed'; this_res <- mets
wb <- openxlsx::createWorkbook()
# creat worksheet
openxlsx::addWorksheet(wb,sprintf('%s', out))
# write data
openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
# create and add a style to the column headers
headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
# style for body
bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
# apply style
addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:(nrow(this_res)+1), cols = 1:ncol(this_res), gridExpand = TRUE)
addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)
# write out
openxlsx::saveWorkbook (wb, file=metabolites_analyzed, overwrite=TRUE)


## finished
print("Done! preprocessing metabolomics data completed.") 
print("Generated excel files with supplementary files in results folder!") 