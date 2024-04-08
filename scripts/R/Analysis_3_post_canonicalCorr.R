#####################################################################
#
# This script is used to conduct post-CCA statistics, including
# 1. general CCA correlation between groups
# 2. Network-wise CCA contribution
# 3. Coupling between CCA map and time delay projection map
#
# Liang Qunjun

library(tidyverse)
library(ggeasy)
library(ggsci)
library(plotly)
library(ggradar)
library(patchwork)

source('scripts/R/function_DrawSurfPlot.R')
source('scripts/R/function_AddNetworkScore.R')
source("scripts/R/function_ColorFroGroup.R")

# path indication
path_cor <-  "inputs/Analysis_CCA_correlation.xlsx"
path_load <- "inputs/Anslysis_CCA_loading.xlsx"
path_td_diff <- "inputs/Analysis2_TDp_difference.xlsx"
net_anna <- readr::read_csv("inputs/Yeo_7net_roiAnnotation.csv", col_names = "net_raw") %>%
  mutate(network = ifelse(str_detect(net_raw, pattern = "Vis"), "VIS",
                          ifelse(str_detect(net_raw, pattern = "SomMot"), "SMN", 
                                 ifelse(str_detect(net_raw, pattern = "DorsAttn"), "DAN",
                                        ifelse(str_detect(net_raw, pattern = "SalVentAttn"), "SAN",
                                               ifelse(str_detect(net_raw, pattern = "Limbic"), "LIM",
                                                      ifelse(str_detect(net_raw, pattern = "Cont"), "FPN", "DMN")))))))
  
########################## visualize CCA correlation ###########################
# load the data
dat_cor <- rio::import(path_cor) 
colnames(dat_cor) <- c("group", colnames(dat_cor)[1:17])
dat_cor <-  dat_cor %>% 
  mutate(subtype = c("LPA","LPA","LPA","LPA","Severity","Severity","Severity")) %>%
  mutate(y = "1")
dat_cor$group <- factor(dat_cor$group, levels = c("saMDD","sdMDD","mdMDD","miMDD","All","severe","modulated"))

# significance of CCA correlation
rho = dat_cor[1,2:18] %>% as.matrix() %>% as.vector()
CCP::p.asym(rho = rho, N = 135, p = 17, q = 400, tstat = "Wilks")

# visualization 
p_dot_cor <- ggplot(dat_cor, aes(x = group, y = y)) +
  geom_point(aes(size = cor01,
                 color = subtype)) +
  geom_text(aes(label = round(cor01, digits = 2))) +
  ggtitle("Correlation for each CCA") +
  scale_size_continuous(range = c(10,25)) +
  bruceR::theme_bruce() + easy_remove_legend() +
  easy_remove_y_axis(what = c("title","ticks","text")) +
  easy_remove_x_axis(what = c("title")) 
p_dot_cor

######################## Statistics for CCA loading ###########################

# load the CCA-Loading
dat_load <- rio::import(path_load) %>% 
  mutate(subtype = c("LPA","LPA","LPA","LPA","Severity","Severity","Severity"), .before = 1) 

# for severe and modulated MDD
p_load_den <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "Severity", group != "All") %>%
  ggplot(aes(x = value, fill = group)) +
  geom_density(alpha = .3, color = "white") +
  xlim(c(-.2,.2)) +
  scale_fill_manual(values=c("steelblue","#F3C846")) +
  theme_classic()
p_load_den

## test the distributions
dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "Severity", group != "All") %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.115, p-value = 0.01008
# alternative hypothesis: two-sided

# for different LPA groups --------------------------------------------------
dat_tmp <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA") 
dat_tmp$group <- factor(dat_tmp$group, levels = c("saMDD","sdMDD","mdMDD","miMDD"))
p_load_den_lpa <- ggplot(dat_tmp, aes(x = value, fill = group)) +
  geom_density(alpha = .5, color = "white") +
  scale_color_manual(values=c(ColorForGroup("saMDD"),ColorForGroup("sdMDD"), 
                              ColorForGroup("miMDD"), ColorForGroup("mdMDD"))) +
  #xlim(c(-.2,.2)) +
  theme_classic()
p_load_den_lpa

## test for the distributions
res_ks1_sasd <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("saMDD", "sdMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.0925, p-value = 0.06526
# alternative hypothesis: two-sided

res_ks1_mimd <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("miMDD", "mdMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.0475, p-value = 0.7576
# alternative hypothesis: two-sided

res_ks1_sdmd <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("sdMDD", "mdMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.1675, p-value = 2.674e-05
# alternative hypothesis: two-sided

res_ks1_sdmi <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("sdMDD", "miMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.17, p-value = 1.908e-05
# alternative hypothesis: two-sided

res_ks1_samd <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("saMDD", "mdMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.1275, p-value = 0.002999
# alternative hypothesis: two-sided

res_ks1_sami <- dat_load %>% pivot_longer(cols = 3:402, names_to = "region") %>%
  filter(subtype == "LPA", group %in% c("saMDD", "miMDD")) %>%
  ks.test(formula = value ~ group, data = .)
# Asymptotic two-sample Kolmogorov-Smirnov test
# data:  value by group
# D = 0.115, p-value = 0.01008
# alternative hypothesis: two-sided

p.adjust(p = c(res_ks1_sasd$p.value, res_ks1_mimd$p.value, 
               res_ks1_sdmd$p.value, res_ks1_sdmi$p.value, 
               res_ks1_samd$p.value, res_ks1_sami$p.value), method = "fdr") %>%
  round(digits = 3)
# 0.078 0.758 0.000 0.000 0.006 0.015

######################### collect in networks ################################
dat_load_abs <- dat_load %>% mutate(across(where(is.numeric), abs))
dat_load_net <- AddNetworkScore_Schaefer(dat_load_abs, net_anna, prefix = 2)
# loading for severe and modulated MDDs
dat_tmp <- dat_load_net %>% filter(subtype == "Severity", group != "All") %>% 
  select(2, starts_with("td_"))
colnames(dat_tmp) <- str_replace(colnames(dat_tmp), pattern = "td_", replacement = "")
p_radar_severe <- ggradar(dat_tmp, grid.mid = 0.03, grid.max = 0.06, values.radar = c("0","0.03","0.06"),
        group.line.width = 1, group.point.size = 3, 
        legend.text.size = 10, legend.title = "subtype")+
  scale_color_manual(values=c("steelblue","#F3C846"))
# loading for LPA groups
dat_tmp <- dat_load_net %>% filter(subtype == "LPA") %>% 
  select(2, starts_with("td_"))
colnames(dat_tmp) <- str_replace(colnames(dat_tmp), pattern = "td_", replacement = "")
p_radar_lpa <- ggradar(dat_tmp, grid.mid = 0.03, grid.max = 0.07, 
        values.radar = c("0","0.03","0.07"),
        group.line.width = 1, group.point.size = 3, 
        legend.text.size = 10, legend.title = "subtype") +
  scale_color_manual(values=c(ColorForGroup("mdMDD"),ColorForGroup("miMDD"), 
                              ColorForGroup("saMDD"), ColorForGroup("sdMDD")))

######################### combine all plots ##################################
layout_use <- "
AAAA
BBCC
BBCC
DDEE
DDEE
"
p_final <- p_dot_cor + p_load_den + p_load_den_lpa + p_radar_severe + p_radar_lpa +
  plot_layout(design = layout_use, guides = 'collect') +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 15, face = "bold"))

ggsave(plot = p_final, filename = "outputs/Fig4.png", width = 10, height = 10)
ggsave(plot = p_final, filename = "outputs/Fig4.tiff", width = 10, height = 10)
