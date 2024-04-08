#########################################################################
#
# Gene expression across different MDD groups
# The expression of Schaefer 400 was obtained from python package abagen
#
#
# Liang Qunjun 2023-12-14

library(tidyverse)
library(ggeasy)
library(pls)

# load the CCA loading map
path_load <- "inputs/Anslysis_CCA_loading.xlsx"
path_gene <- "inputs/Schaefer400_gene_expression.csv"

# load the CCA-Loading
dat_load <- rio::import(path_load) %>% 
  mutate(subtype = c("LPA","LPA","LPA","LPA","Severity","Severity","Severity"), .before = 1) 
colnames(dat_load) <- c("subtype","group", paste0("R",formatC(1:400, width = 3, flag = "0")))

# load the gene expression
dat_gene <- rio::import(path_gene)

set.seed(1234)
n_comp = 6
vali_method = "CV" # for "CV" and "LOO"

############################## for saMDD ###################################
dat_load_tmp <- dat_load %>% filter(group == "saMDD") %>% .[,3:402] %>% as.matrix()
dat_gene_load <- dat_gene %>% mutate("y" = c(dat_load_tmp))

model_psl <- plsr(formula = y ~ ., ncomp = n_comp, 
                  data = dat_gene_load[,2:ncol(dat_gene_load)], 
                  validation = vali_method)
pls::explvar(model_psl)
# Comp 1    Comp 2    Comp 3    Comp 4    Comp 5    Comp 6 
# 10.478221  8.814700  5.923439  6.079632  8.672237  3.799517 
comp_tmp <- model_psl$projection %>% as.data.frame() %>% 
  mutate(gene = rownames(model_psl$projection))
comp_tmp2 <- comp_tmp %>% select(gene, `Comp 1`)
comp_tmp3 <- comp_tmp2[order(comp_tmp2$`Comp 1`, decreasing = T),] %>% select(1)
write_csv(comp_tmp3, col_names = F,
            file = "inputs/Analysis_CCA_gene_rank_saMDD.csv")

############################## for sdMDD ###################################
dat_load_tmp <- dat_load %>% filter(group == "sdMDD") %>% .[,3:402] %>% as.matrix()
dat_gene_load <- dat_gene %>% mutate("y" = c(dat_load_tmp))

model_psl <- plsr(formula = y ~ ., ncomp = n_comp, 
                  data = dat_gene_load[,2:ncol(dat_gene_load)], 
                  validation = vali_method)
pls::explvar(model_psl)
# Comp 1    Comp 2    Comp 3    Comp 4    Comp 5    Comp 6 
# 11.543480  6.355482  7.809301  7.559835  4.901763  5.710318
comp_tmp <- model_psl$projection %>% as.data.frame() %>% 
  mutate(gene = rownames(model_psl$projection))
comp_tmp2 <- comp_tmp %>% select(gene, `Comp 1`)
comp_tmp3 <- comp_tmp2[order(comp_tmp2$`Comp 1`, decreasing = T),] %>% select(1)
write_csv(comp_tmp3, col_names = F,
            file = "inputs/Analysis_CCA_gene_rank_sdMDD.csv")

############################## for mdMDD ###################################
dat_load_tmp <- dat_load %>% filter(group == "mdMDD") %>% .[,3:402] %>% as.matrix()
dat_gene_load <- dat_gene %>% mutate("y" = c(dat_load_tmp))

model_psl <- plsr(formula = y ~ ., ncomp = n_comp, 
                  data = dat_gene_load[,2:ncol(dat_gene_load)], 
                  validation = vali_method)
pls::explvar(model_psl)
# Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6 
# 7.236657 9.564445 5.003427 7.194833 6.250512 5.440751 
comp_tmp <- model_psl$projection %>% as.data.frame() %>% 
  mutate(gene = rownames(model_psl$projection))
comp_tmp2 <- comp_tmp %>% select(gene, `Comp 2`)
comp_tmp3 <- comp_tmp2[order(comp_tmp2$`Comp 2`, decreasing = T),] %>% select(1)
write_csv(comp_tmp3, col_names = F,
            file = "inputs/Analysis_CCA_gene_rank_mdMDD.csv")

############################## for I-MDD ###################################
dat_load_tmp <- dat_load %>% filter(group == "miMDD") %>% .[,3:402] %>% as.matrix()
dat_gene_load <- dat_gene %>% mutate("y" = c(dat_load_tmp))

model_psl <- plsr(formula = y ~ ., ncomp = n_comp, 
                  data = dat_gene_load[,2:ncol(dat_gene_load)], 
                  validation = vali_method)
pls::explvar(model_psl)
# Comp 1   Comp 2   Comp 3   Comp 4   Comp 5   Comp 6 
# 6.041972 8.802681 8.307047 4.807750 6.128572 8.561763
comp_tmp <- model_psl$projection %>% as.data.frame() %>% 
  mutate(gene = rownames(model_psl$projection))
comp_tmp2 <- comp_tmp %>% select(gene, `Comp 2`)
comp_tmp3 <- comp_tmp2[order(comp_tmp2$`Comp 2`, decreasing = T),] %>% select(1)
write_csv(comp_tmp3, col_names = F,
          file = "inputs/Analysis_CCA_gene_rank_miMDD.csv")
