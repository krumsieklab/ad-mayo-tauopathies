# Generates Table 1, Figures 2, 3, from the manuscript
# Figure 1 and Figure 4 were manually created using Adobe illustrator

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath", 'results.makepath')))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries ----
library(maplet) # omics analysis

library(tidyverse) # tidy
library(magrittr) # %<>%

library(openxlsx) # excel
library(readxl) # multi sheet excel

library(gtsummary) # table
library(ggplot2) # plot

library(UpSetR) # upset
library(eulerr) # venn diagram
library(pheatmap) # heatmaps

library(RColorBrewer) # color palettes

source("custom_functions.R") # customized functions

# input colors
annocols <- brewer.pal(n = 9, name = "Pastel1")
# input files
met_pdata <- 'results/tmp_mayo_metabolomics_processed_data.xlsx'
cer_outfile <- 'results/supplementary_table_X_metabolomics_associations_cer.xlsx'
tcx_outfile <- 'results/supplementary_table_X_metabolomics_associations_tcx.xlsx'
# outfile
table_outfile <- 'results/table1.docx'

# load preprocessed data
D <- mt_load_se_xls(file=met_pdata)
# load results 
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

# Table 1 ----
Y <- colData(D) %>% data.frame() %>% 
  filter(!is.na(APOE)) %>%
  mutate(apoe_genotype= case_when(APOE == '44' ~ 2, APOE == '24' ~ 1,
                                  APOE == '34' ~ 1, TRUE ~ 0))%>%
  select(Sex, AgeAtDeath, Diagnosis, BrainRegion, apoe_genotype) %>%
  filter(Diagnosis!='PA') %>% # removing pathologic aging samples
  mutate(AgeAtDeath=case_when(AgeAtDeath=='90+' ~ '90', TRUE~ AgeAtDeath), 
         AgeAtDeath=as.numeric(as.matrix(AgeAtDeath)),
         grp=paste0(BrainRegion, '-', Diagnosis)
  )

# write out
Y %>% tbl_summary(by=grp, include=c('Sex', 'AgeAtDeath', 'apoe_genotype')) %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = table_outfile)

# Figure 2a ----
# distribution of metabolic classes ----
pie_stats <- # row data to get metabolites and their annotations
  D %>% rowData()%>% data.frame() %>% 
  select(name, SUB_PATHWAY, SUPER_PATHWAY) %>% 
  # rename the ones without annotations to other
  mutate(SUPER_PATHWAY = case_when(is.na(SUPER_PATHWAY)~'Other', TRUE~SUPER_PATHWAY)) %>% 
  # format the data frame
  unique() %>% select(SUPER_PATHWAY) %>% table() %>% reshape2::melt()
# name the columns 
names(pie_stats) <- c('Class', 'value')
# compute proportions
pie_stats %<>% mutate(prop = round(100*(value / sum(value)), 2),      
                      cumulative = cumsum(value),
                      midpoint = cumulative - value / 2, 
                      label = paste0(Class, " ", round(value / sum(value) * 100, 1), "%"))

# factor levels need to be the opposite order of the cumulative sum of the values
pie_stats$Class <- factor(pie_stats$Class, levels=rev(pie_stats$Class))

# pie plot ----
pdf(file="results/Figure2a_pie.pdf", width=3, height=3)
print(ggplot(pie_stats, aes(x = 1, weight = value, fill = Class)) +
        geom_bar(width = 1, position = "stack", size=1.3, color='#FDFEFE') +
        coord_polar(theta = "y") +
        geom_text(aes(x = 1.3, y = midpoint, label = label), size=6) +
        theme_void() + theme(text = element_text(size=15))+ 
        scale_fill_manual(values=rev(annocols))+ theme(legend.position ='none'))
dev.off()

# Figure 2 b ----
# four sets to plot
four_sets_res_list <- list(cer_ad=prim_res$CER_AD %>% filter(significant=='Yes') %>% pull(name), 
                           tcx_ad=prim_res$TCX_AD %>% filter(significant=='Yes') %>% pull(name), 
                           cer_psp=prim_res$CER_PSP %>% filter(significant=='Yes') %>% pull(name), 
                           tcx_psp=prim_res$TCX_PSP %>% filter(significant=='Yes') %>% pull(name)
)
# format the data
tmp <- reshape2::melt(lapply(four_sets_res_list, length))
names(tmp) <- c('pheno_num', 'pheno')
tmp$pheno <- factor(tmp$pheno, levels=rev(c('cer_ad', 'cer_psp', 'tcx_ad', 'tcx_psp')))
pdf(file='results/Figure2b_barplot.pdf', height=3.5, width=5)
plot(
  ggplot(tmp, aes(y=pheno_num, x=pheno, fill=pheno, label=pheno_num)) + 
    geom_bar(stat="identity", colour='black') +
    geom_text(size = 5)+
    ylab("#significant associations") + 
    xlab("") + 
    theme_classic() + 
    theme(text = element_text(size=15))+
    ggtitle ("") +
    scale_fill_manual(values=c("#FFFFFF",'#566573', "#FFFFFF", '#566573'))+
    coord_flip(clip = 'off') + theme(legend.position = 'none')
)

dev.off()

# Figure 2 c and d -----
# venns ----
# CER 
pdf(file='results/Figure2c_venn.pdf', height=3, width=3)
print(plot(venn(list(cer_ad=prim_res$CER_AD %>% filter(significant=='Yes') %>% pull(name), 
                     cer_psp=prim_res$CER_PSP %>% filter(significant=='Yes') %>% pull(name))),
           fill=c('#566573', "#FFFFFF", '#CCCCCC'), colour='black', edges=T))
dev.off()
# TCX
pdf(file='results/Figure2d_venn.pdf', height=3, width=3)
print(plot(venn(list(tcx_ad=prim_res$TCX_AD %>% filter(significant=='Yes') %>% pull(name), 
                     tcx_psp=prim_res$TCX_PSP %>% filter(significant=='Yes') %>% pull(name))),
           fill=c('#566573', "#FFFFFF", '#CCCCCC'),edges=T))
dev.off()

# Figure 3 ----
# counts per pathway per group per direction
tmp <- prim_res %>% bind_rows() %>% dplyr::rename(Phenotype=group)
plot_mat <- get_pathway_heat_df(met_stats=tmp, 
                                annocols=annocols,grp1tag='CER_AD', grp2tag='CER_PSP') 

# formatting data for heatmap
plot_mat %<>% select(SUB_PATHWAY, variable, value) %>% 
  # remove the cushions that are remnanats of barplots
  filter(!variable %in%c('tmp', 'co_tmp')) %>% 
  # spread based on sub_pathway
  spread(key=SUB_PATHWAY, value=value) %>% 
  # transform such that pathways are rows
  t() %>% data.frame() %>%
  mutate_all(as.matrix) %>% mutate_all(as.numeric)

# remove first row, contains non info
plot_mat <- plot_mat[-1, ]

# order the data according to the super pathways
tmp_data <- bind_cols(plot_mat, SUB=row.names(plot_mat))
tmp_data <- tmp_data %>%
  left_join(tmp %>% select("SUB_PATHWAY", "SUPER_PATHWAY"), by=c('SUB'='SUB_PATHWAY')) %>% 
  unique()
tmp_data <- tmp_data[order(tmp_data$SUPER_PATHWAY, tmp_data$SUB), ]
annot_mat <- tmp_data %>% select(SUPER_PATHWAY) %>% dplyr::rename(SUP=SUPER_PATHWAY)

# color list of annotations i.e super pathways
annocols_list <- annocols[-length(annocols)] # # defined above pastel colors
names(annocols_list) <- unique(annot_mat$SUP)
annocols_list <- list(SUP=annocols_list)
row_names <- tmp_data$SUB
plot_mat_2 <- tmp_data %>% select(-c(SUPER_PATHWAY, SUB))

# format names of the pathways
row_names <- sub(" Metabolism", "", row_names)
row_names <- sub("and", "&", sub('Fatty Acid', 'FA', row_names))
# order: AD and PSP, positive and then negative
rownames(annot_mat) <- rownames(plot_mat_2) <- row_names
plot_mat_2 <- plot_mat_2 [, c(4,3,2, 1)]

# annotation of unique to AD/PSP or common to both
ad_levels <- rowSums(plot_mat_2[, 1:2])
psp_levels <- rowSums(plot_mat_2[, 3:4])

# add a column in the annotation matrix about the set
annot_mat$set <- 'only AD'
annot_mat$set[intersect(which(psp_levels>0), which(ad_levels==0))] <- 'only PSP'
annot_mat$set[intersect(which(psp_levels>0), which(ad_levels>0))] <- 'AD & PSP'
# add colurs
annocols_list <- list(SUP=annocols_list$SUP, 
                      set=c('#666666','#566573',  '#ffffff'))
names(annocols_list$set) <- c('AD & PSP', 'only AD', 'only PSP')

pdf(height=15, width = 10, file='results/Figure3_heatmap.pdf')
pheatmap(plot_mat_2, cluster_rows = F, cluster_cols = F, 
         display_numbers = T,number_format = '%.0f',
         cellwidth = 10, cellheight = 10, 
         gaps_col = c(2),
         annotation_row = annot_mat, 
         annotation_colors = annocols_list,
         color='white', legend=F)
dev.off()

# Figure 3 venn ----

psp_paths <- tmp %>% filter(significant=='Yes' & Phenotype=='CER_PSP') %>% filter(!is.na(SUB_PATHWAY)) %>% pull(SUB_PATHWAY) %>% unique()
ad_paths <- tmp %>% filter(significant=='Yes' & Phenotype=='CER_AD') %>% filter(!is.na(SUB_PATHWAY)) %>%pull(SUB_PATHWAY) %>% unique()

pdf(file='results/Figure3_venn.pdf', height=3, width=3)
# cer (ad vs ctrl) vs cer (psp vs ctrl) pathways
print(plot(venn(list(cer_ad=ad_paths, 
                     cer_psp=psp_paths)),
           fill=c('#566573', "#FFFFFF", '#CCCCCC'),edges=T))
dev.off()

### finished ------
print("Done!") 
print("Generated word with table 1 and panels of figure 2 and figure 3 in results folder!") 
