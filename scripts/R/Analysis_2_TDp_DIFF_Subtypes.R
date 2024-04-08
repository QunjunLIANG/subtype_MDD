##############################################################################
#
# Time delay projection map for HC contrast 
#
# This script is used to examine the HC-contrast time delay estimation and 
# the diverge pattern of the HC-contrast maps among the MDD subtypes.
#
# The analysis includes:
# 1. t test for HC and each subtype's TDp map.
# 2. F test for network-based difference.
# 3. F test for cohen'd gap among MDD subtypes.
#
# Liang Qunjun 2023/06/07

library(tidyverse)
library(bruceR)
library(ggstatsplot)
library(ggridges)
library(psych)
library(RColorBrewer)
library(emmeans)
library(ggeasy)
library(ggsci)
library(cowplot)
library(purrr)
library(scales)
library(ggsignif)
library(patchwork)
library(sjPlot)
source('scripts/R/function_ObtainNetID.R')
source('scripts/R/function_ObtainBrainData_individual.R')
source('scripts/R/function_AddNetworkScore.R')
source("scripts/R/function_PermuTestBrain.R")
source("scripts/R/function_DrawSurfPlot.R")
source("scripts/R/function_ColorFroGroup.R")

# indicate the path to the inputs and participant information
path_td <- "inputs/timeseries_Schaefer_typical/time_lag_estimation/"
sbj_info <- rio::import('inputs/Analysis1_subject_table.xlsx')
net_anna <- readr::read_csv("inputs/Yeo_7net_roiAnnotation.csv", col_names = "net_raw") %>%
  mutate(network = ifelse(str_detect(net_raw, pattern = "Vis"), "VIS",
                          ifelse(str_detect(net_raw, pattern = "SomMot"), "SMN", 
                                 ifelse(str_detect(net_raw, pattern = "DorsAttn"), "DAN",
                                        ifelse(str_detect(net_raw, pattern = "SalVentAttn"), "SAN",
                                               ifelse(str_detect(net_raw, pattern = "Limbic"), "LIM",
                                                      ifelse(str_detect(net_raw, pattern = "Cont"), "FPN", "DMN")))))))
outfile_name <- "inputs/Analysis2_TDp_collection.xlsx"

if (file.exists(outfile_name)) {
  dat_td <- rio::import(outfile_name)
}else{
  # obtain all file list
  file_names <- list.files(path_td, pattern = 'sub-[0-9]*_projection_map_weighted', full.names = T)
  ## file in used
  file_use <- file_names[grep(paste(sbj_info$participant_id, collapse = "|"), file_names)]
  dat_td_raw <- map_dfr(data.frame(file_use), ObtainBrainData_individual, .progress = T)
  colnames(dat_td_raw)[2:ncol(dat_td_raw)] <- paste0("R",formatC(1:400, width = 3, flag = "0"))
  ## add network-level mean TDp
  dat_td_net <- AddNetworkScore_Schaefer(dat_td_raw, net_anna)
  ## add group index
  dat_td <- dat_td_net %>%
    left_join(sbj_info)
}
## factorize the LPAgroup
dat_use <- dat_td

############################################################
#
# Construct HC-contrast maps
set.seed(1234)

res_samdd <-  TtestBrain(dat_use, subtype1 = "saMDD")
res_sdmdd <-  TtestBrain(dat_use, subtype1 = "sdMDD")
res_mdmdd <-  TtestBrain(dat_use, subtype1 = "mdMDD")
res_mimdd <-  TtestBrain(dat_use, subtype1 = "miMDD")

###########################################################
#
# visualize the t statistic on brain surface

d_collect <- data.frame(samdd = res_samdd$d_value, sdmdd = res_sdmdd$d_value,
                        mdmdd = res_mdmdd$d_value, mimdd = res_mimdd$d_value)

name_collect <- c("saMDD","sdMDD","mdMDD","miMDD")
p_tdiff <- DrawMultipleSurfaceOnSchaefer400(dat_value = d_collect, 
                                 dat_name = name_collect, p_title = "Time shifts in subtypes",
                                 legend_title = "Cohen's d", legent_pos = "left")

ggsave(plot = p_tdiff, filename = "outputs/Anay2_d_difference.png", width = 9, height = 7)

##############################################################
#
# Network-based difference in HC contrasts
#
## test TD difference between saMDD and HC
dat_use %>% filter(subtype_lpa == "saMDD" | subtype_lpa == "HC") %>%
  select(participant_id, subtype_lpa, starts_with("td_"), age, gender) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "subtype_lpa", within = 'network',
                 covariate = c("age","gender"))%>%
  bruceR::EMMEANS(effect = "subtype_lpa", by = "network", p.adjust = "fdr")
# ───────────────────────────────────────────────────────────────────────────────────
# Contrast "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────────────
# saMDD - HC    td_DAN    0.015 (0.010) 118  1.444  .151      0.190 [-0.070,  0.450]
# saMDD - HC    td_DMN    0.001 (0.008) 118  0.173  .863      0.018 [-0.186,  0.222]
# saMDD - HC    td_FPN    0.007 (0.011) 118  0.648  .518      0.091 [-0.188,  0.370]
# saMDD - HC    td_LIM   -0.005 (0.021) 118 -0.262  .794     -0.070 [-0.597,  0.457]
# saMDD - HC    td_SAN    0.012 (0.011) 118  1.091  .277      0.161 [-0.131,  0.453]
# saMDD - HC    td_SMN   -0.018 (0.008) 118 -2.199  .030 *   -0.235 [-0.446, -0.023]
# saMDD - HC    td_VIS   -0.009 (0.010) 118 -0.862  .390     -0.114 [-0.376,  0.148]
# ───────────────────────────────────────────────────────────────────────────────────

## test TD difference between miMDD and HC
dat_use %>% filter(subtype_lpa == "miMDD" | subtype_lpa == "HC") %>%
  select(participant_id, subtype_lpa, starts_with("td_"), age, gender) %>%
  pivot_longer(cols = 3:9, names_to = "network", values_to = "td") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td", 
                 between = "subtype_lpa", within = 'network',
                 covariate = c("age","gender"))%>%
  bruceR::EMMEANS(effect = "subtype_lpa", by = "network", p.adjust = "fdr")
# ───────────────────────────────────────────────────────────────────────────────────
# Contrast "network" Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────────────
# miMDD - HC    td_DAN    0.033 (0.011) 117  3.016  .003 **    0.434 [ 0.149, 0.719]
# miMDD - HC    td_DMN   -0.002 (0.008) 117 -0.278  .782      -0.029 [-0.232, 0.175]
# miMDD - HC    td_FPN   -0.001 (0.011) 117 -0.111  .911      -0.016 [-0.309, 0.276]
# miMDD - HC    td_LIM   -0.037 (0.020) 117 -1.861  .065 .    -0.490 [-1.011, 0.031]
# miMDD - HC    td_SAN   -0.016 (0.012) 117 -1.353  .179      -0.207 [-0.510, 0.096]
# miMDD - HC    td_SMN   -0.004 (0.008) 117 -0.482  .631      -0.053 [-0.273, 0.166]
# miMDD - HC    td_VIS   -0.002 (0.010) 117 -0.200  .842      -0.025 [-0.275, 0.225]
# ───────────────────────────────────────────────────────────────────────────────────

##############################################################
#
# Network-based difference in TDp - visualization
#

p_diff_samdd <- dat_use %>% filter(subtype_lpa == "saMDD" | subtype_lpa == "HC") %>%
  select(participant_id, subtype_lpa, td_SMN) %>%
  ggplot(aes(x = subtype_lpa, y = td_SMN, fill = subtype_lpa)) +
  geom_violin(width = .5, alpha = .5, color = "white") +
  geom_boxplot(width = .25, alpha = .7) +
  coord_flip() +
  ylab("Time delay estimation") + ggtitle("Somatomotor netowrk") +
  scale_fill_manual(values = c("steelblue",ColorForGroup("saMDD"))) +
  theme_classic() + easy_remove_legend() + easy_remove_y_axis(what = "title") +
  easy_text_size(12)

p_diff_mimdd <- dat_use %>% filter(subtype_lpa == "miMDD" | subtype_lpa == "HC") %>%
  select(participant_id, subtype_lpa, td_DAN) %>%
  ggplot(aes(x = subtype_lpa, y = td_DAN, fill = subtype_lpa)) +
  geom_violin(width = .5, alpha = .5, color = "white") +
  geom_boxplot(width = .25, alpha = .7) +
  coord_flip() +
  ylab("Time delay estimation") + ggtitle("Dorsal attention netowrk") +
  scale_fill_manual(values = c("steelblue",ColorForGroup("saMDD"))) +
  theme_classic() + easy_remove_legend() + easy_remove_y_axis(what = "title") +
  easy_text_size(12)

##############################################################
#
# HC-contrast map difference, evaluated by Cohen's d distribution
#

# obtain effect size - conhens'd
dat_eff <- data.frame(saMDD = res_samdd$d_value, sdMDD = res_sdmdd$d_value,
                       miMDD = res_mimdd$d_value, mdMDD = res_mdmdd$d_value)
dat_eff_long <- dat_eff %>% pivot_longer(cols = 1:4, names_to = "group")
dat_eff_long$group <- factor(dat_eff_long$group, levels = c("saMDD","sdMDD","mdMDD","miMDD"))

dat_eff_long %>% mutate(value_abs = abs(value)) %>%
  bruceR::MANOVA(dv = "value_abs", between = "group") %>%
  bruceR::EMMEANS(effect = "group", p.adjust = 'fdr')
# ─────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.   df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────
# sdMDD - saMDD   -0.045 (0.009) 1596 -4.955 <.001 *** -0.350 [-0.537, -0.164]
# mdMDD - saMDD   -0.023 (0.009) 1596 -2.502  .021 *   -0.177 [-0.364,  0.010]
# mdMDD - sdMDD    0.022 (0.009) 1596  2.453  .021 *    0.173 [-0.013,  0.360]
# miMDD - saMDD   -0.004 (0.009) 1596 -0.467  .640     -0.033 [-0.220,  0.154]
# miMDD - sdMDD    0.041 (0.009) 1596  4.488 <.001 ***  0.317 [ 0.131,  0.504]
# miMDD - mdMDD    0.019 (0.009) 1596  2.035  .050 .    0.144 [-0.043,  0.331]
# ─────────────────────────────────────────────────────────────────────────────

p_eff_diff <- ggplot(dat_eff_long, aes(y = abs(value), x = group)) +
  geom_violin(width = .79) +
  geom_point(position= position_jitter(.11), 
             aes(color = group), size = 2.8, alpha = .2) +
  geom_boxplot(width = .4, color = "black", 
               outlier.alpha = 0, fill = NA) +
  geom_signif(annotations = c("***","*","***","*","*"),
              y_position = c(.7,.72,.74,.76,.78), textsize = 6,vjust = .8,
              xmin = c(1,1,2,2,3), xmax = c(2,3,3,4,4), tip_length = 0) +
  scale_color_manual(values=c(ColorForGroup("saMDD"),ColorForGroup("sdMDD"), 
                              ColorForGroup("miMDD"), ColorForGroup("mdMDD"))) +
  ylab("Cohen's d") + 
  theme_classic() + easy_remove_x_axis(what = "title") + easy_text_size(13) +
  easy_add_legend_title("Subtype")
p_eff_diff

##############################################################
#
# Combine the plots
#

layout_use <- "
AAAA
AAAA
AAAA
BBDD
CCDD
"
p_final <- p_tdiff +  p_diff_samdd + p_diff_mimdd + p_eff_diff +
  patchwork::plot_annotation(tag_levels = "A") + 
  patchwork::plot_layout(design = layout_use, guides = "collect") &
  theme(plot.tag = element_text(size = 17, face = "bold"))
p_final
ggsave(plot = p_final, filename = "outputs/Fig2.png", width = 10, height = 11)
ggsave(plot = p_final, filename = "outputs/Fig2.tiff", width = 10, height = 11)

##############################################################
#
# Export the new subject information
#

rio::export(dat_use, file = outfile_name)

dat_d_collect <- d_collect
colnames(dat_d_collect) <- c("saMDD","sdMDD","mdMDD","miMDD")
rio::export(dat_d_collect, file = "inputs/Analysis2_TDp_difference.xlsx")
