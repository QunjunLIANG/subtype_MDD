###########################
#
# This function is used to conducted 
# permutation test for group difference 
# time delay estimation between HC
# and subtype, or between two subtypes.
#
# R packages coin and tidyverse are essential 
# for this function.
#
# Liang Qunjun 

TtestBrain <- function(data, subtype1, subtype2="HC"){
  set.seed(1234)
  library(tidyverse)
  # filter the data
  dat_test <- data %>% filter(subtype_lpa == subtype2 | subtype_lpa == subtype1)
  dat_test$subtype_lpa <- factor(dat_test$subtype_lpa)
  
  regions <- c()
  pvalues <- c()
  tvalues <- c()
  dvalues <- c()
  for (region in 1:400) {
    dat_tmp <- dat_test %>% select(1 + region, subtype_lpa)
    region_tmp <- colnames(dat_tmp)[1] # obtain this region
    
    colnames(dat_tmp)[1] <- "region" # modify the colname
    test <- bruceR::TTEST("region", x = "subtype_lpa", data = dat_tmp)
    # obtain the statistics
    p_tmp <- test$pval
    t_tmp <- test$t
    d_tmp <- test$Cohen_d
    # add to the sequence
    regions <- c(regions, region_tmp)
    pvalues <- c(pvalues,p_tmp)
    tvalues <- c(tvalues,t_tmp)
    dvalues <- c(dvalues, d_tmp)
  }
  dat_out <- data.frame(region = regions, t_value = tvalues, p_value = pvalues, d_value = dvalues)
  
  return(dat_out)
}

PermuTestBrain <- function(data, subtype1, subtype2="HC"){
  set.seed(1234)
  library(tidyverse)
  library(coin)
  # filter the data
  dat_test <- data %>% filter(subtype_lpa == subtype2 | subtype_lpa == subtype1)
  dat_test$subtype_lpa <- factor(dat_test$subtype_lpa)
  
  regions <- c()
  pvalues <- c()
  tvalues <- c()
  for (region in 1:400) {
    dat_tmp <- dat_test %>% select(1 + region, subtype_lpa)
    region_tmp <- colnames(dat_tmp)[1] # obtain this region
    
    colnames(dat_tmp)[1] <- "region" # modify the colname
    test <- coin::oneway_test(region ~ subtype_lpa, # do the test
                        data = dat_tmp, distribution = approximate(nresample = 5000))
    # obtain the statistics
    p_tmp <- coin::pvalue(test)[1] 
    t_tmp <- coin::statistic(test)
    # add to the sequence
    regions <- c(regions, region_tmp)
    pvalues <- c(pvalues,p_tmp)
    tvalues <- c(tvalues,t_tmp)
  }
  dat_out <- data.frame(region = regions, t_value = tvalues, p_value = pvalues)
  
  return(dat_out)
}



