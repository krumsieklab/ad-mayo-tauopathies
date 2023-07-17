# Script to calculate ration of metabolites from polyamine metabolite associated with diagnosis
# for two brain regions - tcx and cer and two diagnosis - AD and PSP
# Generates supplementary tables with statistical results of the association analysis

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath")))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# selected pathwayas
selected_pathways <- c("Polyamine Metabolism")
# input files
met_pdata <- 'results/tmp_mayo_metabolomics_processed_data.xlsx'

# output files
outfile <- 'results/supplementary_table_X_ratio_associations.xlsx'

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
# filter for selected pathways
D1 <- D1 %>% mt_modify_filter_features(filter = SUB_PATHWAY%in%selected_pathways)

# compute ratios
# assay data
X <- D1 %>% assay() %>% data.frame()
# unlog2
X <- 2^X
# create ratios
tmp <- list(); k=1; namen <- NULL
# loop over metabolites
for(i in c(1:nrow(X))){
  # skip the last metabolite
  if(i<nrow(X)){
    # second loop over metabolites for pairwise ratios  
    for(j in (c((i+1):nrow(X)))){
      # ratio and then log2
      tmp[[k]] <- log2((X[i, ]) / (X[j, ]))
      # increase index
      k <- k+1
      # save a new name for the ratio
      namen <- c(namen, paste(rowData(D1)$name[i], rowData(D1)$name[j], sep="_"))
    }
  }
}
# data frame with ratios
ratios <- do.call(rbind, tmp) %>% data.frame()
rownames(ratios) <- paste0("ratio", 1:nrow(ratios))
colnames(ratios) <- colnames(D1)
# new summarized experiment with ratios
D1 <- SummarizedExperiment(assay=ratios, rowData=data.frame(name=rownames(ratios), fullname=namen), colData=colData(D1))

res_list <- list()# empty output list
subsets <- c('AD', 'PSP')# tauopathies
conf_form <- 'Sex + AgeAtDeath + apoe_genotype'
# loop over brain regions
for (brain_region in c('CER')){
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
# workbook
wb <- openxlsx::createWorkbook()
#loop over brain regions
for(brainregion in names(res_list)){

  #loop over per diagnosis
  for(diagnosis in names(res_list[[brainregion]])){
    out <- sprintf('%s_%s', brainregion, diagnosis)
    # create worksheet
    openxlsx::addWorksheet(wb,sprintf('%s', out))
    this_res <- res_list[[brainregion]][[diagnosis]] %>% .[order(.$significant, decreasing = T), ] %>%
      select(-analyte) %>% select(fullname, outcome, covariates, significant, p_value, adj_p1, everything())
    
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

}
# outfile
openxlsx::saveWorkbook (wb, file=outfile, overwrite=TRUE)  
## finished
print("Done! ratio analysis for both taupathies and brain regions.") 
print("Generated excel file with supplementary table in results folder!") 
