# Customized functions used in other scripts.

# This script includes the following functions which can
# be used for two stage analysis
# (a) mayo_analysis
# for any association analysis 
# (b) association_analysis
# for 

# two stage association analysis
mayo_analysis <- function(D, outcome = 'Diagnosis', 
                          outcome_type =  'twofactor', 
                          conf_formula = NULL,padjcut,  pcut, D0, stat_name){
  
  #model 1
  res_m1 <- association_analysis(D, outcome = outcome, 
                                 outcome_type =  outcome_type)
  # model 2
  res_m2 <- association_analysis(D, outcome = outcome, 
                                 outcome_type = outcome_type, 
                                 conf_formula = conf_formula)
  # significant mets
  sigs <- intersect(res_m1 %>% filter(adj_p<=padjcut) %>% pull(name), 
                    res_m2 %>% filter(p_value<=pcut) %>% pull(name))
  # final result table
  res <- left_join(res_m2, res_m1 %>% 
                     select(name, estimate, std_error, statistic, p_value, adj_p) %>% 
                     dplyr::rename(estimate1=estimate, std_error1=std_error,
                                   statistic1=statistic, p_value1=p_value, adj_p1=adj_p), by='name') %>% 
    mutate(significant= case_when(name%in%sigs ~ 'Yes', TRUE ~ 'No'))
  # return res table
  return(res)
}
# statistical associations
association_analysis <- function(D,  # SE
                                 outcome, # outcome to be used for response 
                                 outcome_type, # outcome type - numeric/binary/ordinal
                                 conf_formula=NULL, # confounders to be corrected for
                                 int_w_analyte=NULL, all_vals=F){ # if there is an interaction term
  
  # merge data with sample info
  Ds <- D %>% maplet:::mti_format_se_samplewise() # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  # turn the outcome variable into factor
  Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
  # metabolites in data
  mets <- D %>% assay() %>% rownames()
  # loop over metabolites
  univ_stats <- lapply(mets, function(x){
    
    # formula for this metabolite
    if(is.null(int_w_analyte)==F){
      this_formula <- as.formula(glue('{outcome} ~ {x} + {conf_formula} + {int_w_analyte}*{x}'))
    } else if (is.null(conf_formula)==F){
      this_formula <- as.formula(glue('{outcome} ~ {x} + {conf_formula}'))  
    } else{
      this_formula <- as.formula(glue('{outcome} ~ {x}'))  
    }
    
    # univariate analysis of numeric outcomes
    if(outcome_type=='numeric'){
      # turn this outcome variable into numeric 
      Ds %<>% mutate(!!sym(outcome) := as.numeric(as.matrix(!!sym(outcome))))
      # linear regression
      this_fit <- lm(this_formula, data=Ds)
      # univariate tests of binary outcomes
    } else if (outcome_type=='twofactor'){
      # turn this outcome variable into factor 
      Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
      
      # logistic regression
      this_fit <- glm(this_formula, data=Ds, family=binomial(link='logit'))
    }
    # results summary
    this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
    names(this_res) <- c('estimate', 'std_error', 'statistic', 'p_value')
    this_res <- this_res [2, ]
    this_res %<>% mutate(analyte=rownames(this_res),
                         outcome=outcome, 
                         covariates=conf_formula)
    # order output columns
    this_res %<>% select(analyte, outcome, everything())
    return(this_res)
  }) %>% # create data from of results
    do.call(rbind,.) %>% data.frame() %>% 
    # bind rowData
    bind_cols(D %>% rowData() %>% data.frame() %>% select(name, everything()))
  
  # adjust pvalues
  univ_stats %<>% mutate(adj_p = p.adjust(p_value, method='BH'))
  # order by adjusted p-values
  univ_stats <- univ_stats[order(univ_stats$adj_p), ]
  # return results
  return(univ_stats)
}

# counts per pathway per group per direction
get_pathway_heat_df <- function(met_stats, annocols, grp1tag, grp2tag){
  
  grp1 <- met_stats %>% filter(significant=='Yes' & Phenotype==grp1tag) %>%
    mutate(SUB_PATHWAY=case_when(is.na(SUB_PATHWAY)  ~ 'Other', TRUE~SUB_PATHWAY))
  grp1 <- dplyr::left_join(grp1, grp1 %>% dplyr::group_by(SUB_PATHWAY) %>%
                             dplyr::count(name) %>% dplyr::group_by(SUB_PATHWAY) %>% 
                             dplyr::count(SUB_PATHWAY), by='SUB_PATHWAY') %>% 
    dplyr::rename(path_wt_grp1=n) %>%
    dplyr::left_join(grp1 %>% group_by(SUB_PATHWAY) %>% 
                       count(statistic<0) %>% 
                       filter(`statistic < 0` ==TRUE) %>% 
                       select(- `statistic < 0`), by='SUB_PATHWAY') %>% 
    dplyr::rename(low_grp1=n) %>%
    dplyr::left_join(grp1 %>% group_by(SUB_PATHWAY) %>% 
                       count(statistic>0) %>% 
                       filter(`statistic > 0` ==TRUE) %>% 
                       select(- `statistic > 0`), by='SUB_PATHWAY') %>% 
    dplyr::rename(high_grp1=n)
  
  grp1 %<>% select(SUB_PATHWAY,path_wt_grp1, low_grp1, high_grp1) %>% unique()
  # metabolomics
  grp2 <- met_stats %>% filter(significant=='Yes' & Phenotype==grp2tag) %>%
    mutate(SUB_PATHWAY=case_when(is.na(SUB_PATHWAY)  ~ 'Other', TRUE~SUB_PATHWAY))
  # format the data
  grp2 <- dplyr::left_join(grp2, grp2 %>% dplyr::group_by(SUB_PATHWAY) %>%
                             dplyr::count(name) %>% dplyr::group_by(SUB_PATHWAY) %>% 
                             dplyr::count(SUB_PATHWAY), by='SUB_PATHWAY') %>% 
    dplyr::rename(path_wt_grp2=n) %>%
    dplyr::left_join(grp2 %>% group_by(SUB_PATHWAY) %>% 
                       count(statistic<0) %>% 
                       filter(`statistic < 0` ==TRUE) %>% 
                       select(- `statistic < 0`), by='SUB_PATHWAY') %>% 
    dplyr::rename(low_grp2=n) %>%
    dplyr::left_join(grp2 %>% group_by(SUB_PATHWAY) %>% 
                       count(statistic>0) %>% 
                       filter(`statistic > 0` ==TRUE) %>% 
                       select(- `statistic > 0`), by='SUB_PATHWAY') %>% 
    dplyr::rename(high_grp2=n)
  
  grp2 %<>% select(SUB_PATHWAY,path_wt_grp2, low_grp2, high_grp2) %>% unique()
  
  tmp <- full_join(grp1, grp2, by='SUB_PATHWAY') %>% unique() %>%
    # remove(other) 
    filter(SUB_PATHWAY!='Other')
  plot_mat <- reshape2::melt(tmp, id=c('SUB_PATHWAY', 'path_wt_grp1', 'path_wt_grp2')) %>%
    mutate(value=as.numeric(as.matrix(value)),
           value=case_when(is.na(value)~0, TRUE~value)) %>% unique()
  plot_mat$grp <- sub('low_', '', sub('high_', '', plot_mat$variable))
  
  y_lab <- '# significant metabolites'
  # adding unique labels so they each group gets its own bar order is important!!!
  test <- bind_rows(plot_mat %>% filter(grp=='grp1') %>% 
                      mutate(label=paste0(SUB_PATHWAY, '_1')),
                    plot_mat %>% filter(grp=='grp2') %>% 
                      mutate(label=paste0(SUB_PATHWAY, '_2')),
                    plot_mat %>% filter(grp=='grp1') %>% 
                      mutate(label=paste0(SUB_PATHWAY, '_3'), 
                             value=0, grp='empty_bars', 
                             variable=case_when(variable=='low_grp1' ~ 'tmp', TRUE~'co_tmp')))
  # level manipulation of labels
  test$SUB_PATHWAY <- factor(test$SUB_PATHWAY, levels=c(rev(levels(factor(test$SUB_PATHWAY)))) )
  # order of groups
  test$variable <- factor(test$variable, levels=rev(c('co_tmp','tmp', 'high_grp1', 'low_grp1', 'high_grp2', 'low_grp2')))
  # alternate sequence that takes two bys
  seq_alt <- function(from, to, by) {
    seq <- cumsum(c(from,rep(by,ceiling((to-from)/sum(by)))))
    return(seq[! seq > to])
  }
  # add x axis values to place the bars in desired distance
  tmp <-  test[order(test$SUB_PATHWAY), ] %>% 
    mutate(idx=sort(rep(seq(0.5, 5000, 0.5)[-seq_alt(3, 10000, c(3, 2))], 2))[1:nrow(test)])
  # xaxis label maniupulation
  pns <- tmp %>% pull(SUB_PATHWAY) %>% sort() %>% # take the pathway names
    sub('..*_1', '', .) %>% # remove the BS labels
    sub('..*_3', '', .) %>% # remove the EB labels
    sub("_2", '', .) %>%  # remove suffix from NS labels
    sub(' pathway', '', .) %>% # remove pathway as suffix
    sub(' signaling', '', .) %>% # remove signaling as suffix
    sub(' Metabolism', '', .) # remove metabolism as suffix
  pns[which(duplicated(pns))] <- ''
  
  return(tmp)
} 

