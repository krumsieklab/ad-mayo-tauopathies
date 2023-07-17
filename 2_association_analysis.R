# Script to calculate metabolites associated with diagnosis
# for two brain regions - tcx and cer and two diagnosis - AD and PSP
# Generates supplementary tables with statistical results of the association analysis

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath")))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# input files
met_pdata <- 'results/tmp_mayo_metabolomics_processed_data.xlsx'

# output files
cer_outfile <- 'results/supplementary_table_X_metabolomics_associations_cer.xlsx'
tcx_outfile <- 'results/supplementary_table_X_metabolomics_associations_tcx.xlsx'

# input parameters
padjcut <- 0.25 # stage 1 adjusted p-value cutoff
pcut <- 0.05 # stage 2 nominal p-value cutoff 

# libraries
library(maplet) # omics analysis 

library(tidyverse) # tidy
library(magrittr) # %<>%

library(glue) # formula

library(openxlsx) # excel

source("custom_functions.R") # customized functions

# result directory - creates if does not exist
dir.create("results/", showWarnings = F, recursive = T)


D <- mt_load_se_xls(file=met_pdata) 
# clean metadata
D1 <- D%>%
  # remove samples without sex, apoe and age info
  mt_modify_filter_samples(!is.na(Sex)) %>% 
  mt_modify_filter_samples(!is.na(AgeAtDeath)) %>% 
  mt_modify_filter_samples(!is.na(APOE)) %>% 
  # count the number of apoe4s and make an ordinal variable
  mt_anno_mutate(anno_type = 'samples', 
                 col_name = 'apoe_genotype',
                 term = case_when(APOE == '44' ~ 2, APOE == '24' ~ 1,
                                  APOE == '34' ~ 1, TRUE ~ 0))%>%
  # numeric sex
  mt_anno_mutate(anno_type = 'samples', 
                 col_name = 'Sex',
                 term = case_when(Sex == 'F' ~ 0, TRUE ~ 1))%>%
  # convert 90+ to 90
  mt_anno_mutate(anno_type = "samples", col_name = 'AgeAtDeath', 
                 term = case_when(AgeAtDeath=='90+' ~ '90', TRUE ~ AgeAtDeath)) %>%
  # convert sex to factors and age to numeric
  mt_anno_mutate(anno_type = "samples", col_name = 'Sex', 
                 term = as.factor(as.factor(Sex))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = 'AgeAtDeath', 
                 term = as.numeric(as.matrix(AgeAtDeath))) %>% 
  mt_anno_mutate(anno_type = "samples", col_name = 'apoe_genotype', 
                 term = as.factor(as.matrix(apoe_genotype))) %>%
  {.}

# compute and save metabolic associations ----
# new SE with reduced rowData
D1 <- SummarizedExperiment::SummarizedExperiment(assays = assay(D1), 
                                                 rowData = (rowData(D1) %>% data.frame() %>%
                                                   select(name, SUB_PATHWAY, SUPER_PATHWAY, HMDb, PUBCHEM, CAS, KEGG, COMP_ID) %>% dplyr::rename(HMDB=HMDb)),
                                                 colData= colData(D1) %>% data.frame(),
                                               metadata = list())

res_list <- list()# empty output list
subsets <- c('AD', 'PSP')# tauopathies
conf_form <- 'Sex + AgeAtDeath + apoe_genotype'
# loop over brain regions
for (brain_region in c('TCX', 'CER')){
  # loop over diagnosis
  for(subset in subsets){
    D2 <- D1 %>%  # select only one brain region
      mt_modify_filter_samples(filter = BrainRegion==brain_region) %>%
      # select samples of a specific diagnosis and controls
      mt_modify_filter_samples(filter = Diagnosis %in% c('Control', subset)) %>%
      # numeric diagnosis
      mt_anno_mutate(anno_type = 'samples', 
                     col_name = 'Diagnosis',
                     term = case_when(Diagnosis == 'Control' ~ 0, TRUE ~ 1)) %>%
      # convert diagnosis to factor
      mt_anno_mutate(anno_type = "samples", col_name = 'Diagnosis', 
                     term = as.factor(as.factor(Diagnosis)))
    # all sample analysis - correct for sex, age, apoe
    # 2 stage analysis
    res <- mayo_analysis (D=D2, outcome = 'Diagnosis', 
                          outcome_type =  'twofactor', 
                          conf_formula = conf_form, padjcut=padjcut, pcut=pcut)
    res_list[[brain_region]][[subset]] <- res
  }
}

# write out

#loop over brain regions
for(brainregion in names(res_list)){
  # workbook
  wb <- openxlsx::createWorkbook()
  #loop over per diagnosis
  for(diagnosis in names(res_list[[brainregion]])){
    out <- sprintf('%s_%s', brainregion, diagnosis)
    # create worksheet
    openxlsx::addWorksheet(wb,sprintf('%s', out))
    this_res <- res_list[[brainregion]][[diagnosis]] %>% .[order(.$significant, decreasing = T), ] %>%
      select(-analyte) %>% select(name, outcome, covariates, significant, p_value, adj_p1, everything())
    
    # write data
    openxlsx::writeData(wb, sprintf('%s', out), this_res, rowNames = F, colNames = T)
    # create and add a style to the column headers
    headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
    # style for body
    bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
    # apply style
    addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:(nrow(this_res)+1), cols = 1:ncol(this_res), gridExpand = TRUE)
    addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)
  }
  # outfile name as per the brain region
  if(grepl('CER', brainregion, ignore.case = T)){
    openxlsx::saveWorkbook (wb, file=cer_outfile, overwrite=TRUE)  
  } else{
    openxlsx::saveWorkbook (wb, file=tcx_outfile, overwrite=TRUE)  
  }
}

## finished
print("Done! metabolomics analysis for both taupathies and brain regions.") 
print("Generated excel file with supplementary tables in results folder!") 
