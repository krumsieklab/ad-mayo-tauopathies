# Creates one supplementary file containing all supplementary tables

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath", 'results.makepath')))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries ----

library(tidyverse) # tidy
library(magrittr) # %<>%

library(openxlsx) # excel
library(readxl) # multi sheet excel

padjcut <- 0.25 # stage 1 adjusted p-value cutoff
# input files
mets <- 'results/supplementary_table_X_metabolites_analyzed.xlsx'
cer_outfile <- 'results/supplementary_table_X_metabolomics_associations_cer.xlsx'
tcx_outfile <- 'results/supplementary_table_X_metabolomics_associations_tcx.xlsx'
ratios <- 'results/supplementary_table_X_ratio_associations.xlsx'

# output file
outfile <- 'results/stable.xlsx'
# load statistical results from analysis of metabolites
prim_res_1 <- cer_outfile %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  map(readxl::read_excel, path = cer_outfile) %>%
  purrr::imap(~mutate(.x, group = .y))
prim_res_2 <- tcx_outfile %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  map(readxl::read_excel, path = tcx_outfile) %>%
  purrr::imap(~mutate(.x, group = .y))
# join cer and tcx results
prim_res <- c(prim_res_1, prim_res_2)

# format the above tables
prim_res <- lapply(prim_res, FUN=function(x) {
  x %<>% select(name, outcome, covariates, significant, "estimate1", "std_error1",
                "statistic1", "p_value1", adj_p1, "estimate", "std_error", "statistic", 
                "p_value", everything()) %>%
    select(-adj_p, -group) %>%
    dplyr::rename(estimate_model1=estimate1, 
                  std_error_model1=std_error1, 
                  statistic_model1=statistic1,
                  p_value_model1=p_value1,
                  adj_p_model1=adj_p1,
                  estimate_model2=estimate, 
                  std_error_model2=std_error, 
                  statistic_model2=statistic,
                  p_value_model2=p_value) %>%
    mutate(estimate_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~estimate_model2), 
           std_error_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~std_error_model2),
           statistic_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~statistic_model2),
           p_value_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~p_value_model2))
  x <- x [order(x$p_value_model2, decreasing = F), ]
  return(x)
}) 
# load statistical results from analysis of ratios
ratio_res <- ratios %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  map(readxl::read_excel, path = ratios) %>%
  purrr::imap(~mutate(.x, group = .y))

# format the above tables
ratio_res <- lapply(ratio_res, FUN=function(x) {
  x %<>% select(fullname, outcome, covariates, significant, "estimate1", "std_error1",
                "statistic1", "p_value1", adj_p1, "estimate", "std_error", "statistic", 
                "p_value", everything()) %>%
    select(-adj_p, -group, -name) %>%
    dplyr::rename(name=fullname, 
                  estimate_model1=estimate1, 
                  std_error_model1=std_error1, 
                  statistic_model1=statistic1,
                  p_value_model1=p_value1,
                  adj_p_model1=adj_p1,
                  estimate_model2=estimate, 
                  std_error_model2=std_error, 
                  statistic_model2=statistic,
                  p_value_model2=p_value) %>%
    mutate(estimate_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~estimate_model2), 
           std_error_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~std_error_model2),
           statistic_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~statistic_model2),
           p_value_model2=case_when(adj_p_model1>padjcut ~ NA_real_, TRUE~p_value_model2))
  x <- x [order(x$p_value_model2, decreasing = F), ]
  return(x)
}) 

# one table with all the tables 
mets  <- mets %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  map(readxl::read_excel, path = mets) %>%
  purrr::imap(~mutate(.x, group = .y))
all_tabs <-NULL
all_tabs <- c(mets, prim_res, ratio_res)
names(all_tabs) <- paste0("stable_", 1:length(all_tabs))
# write out
# workbook
wb <- openxlsx::createWorkbook()
#loop over brain regions
for(i in c(1:length(all_tabs))){
  out <- names(all_tabs)[i]
  # create worksheet
  openxlsx::addWorksheet(wb,sprintf('%s', out))
  this_res <- all_tabs[[i]]
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
# outfile
openxlsx::saveWorkbook (wb, file=outfile, overwrite=TRUE)  

### finished ------
print("Done!") 
print("Generated supplementary file in results folder!") 
