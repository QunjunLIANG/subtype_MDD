##############################################################################
#
# Canonical correlation analysis
#
# This script is used to obtain the correlation between time delay using
# canonical correlation analysis.
#
# We obtained the correlation value and the CCA loading of each region.
#
# Note. The CCA used here is not the traditional CCA, but a regular CCA.
# Please refer to the R package CCA for more details.
# 
# Running of the R. CCA is a costly of time. Therefore, the statistics were
# proceed in the next script. Here, we just obtain the raw-data of R. CCA.
#
# Liang Qunjun 2023-12-11

library(tidyverse)
library(bruceR)
library(ggridges)
library(psych)
library(ggeasy)
library(ggsci)
library(purrr)
library(scales)
library(patchwork)
source('scripts/R/function_DrawSurfPlot.R')
source('scripts/R/function_ObtainNetID.R')
source('scripts/R/function_ObtainBrainData_individual.R')
source('scripts/R/function_AddNetworkScore.R')

# indicate the path to the inputs and participant information
path_td <- "inputs/timeseries_Schaefer_typical/time_lag_estimation/"
path_td_diff <- "inputs/Analysis2_TDp_difference.xlsx"
sbj_info <- rio::import('inputs/Analysis1_subject_table.xlsx')
net_anna <- rio::import("inputs/Yeo_net_identity_7.csv") %>% 
  mutate(network = ifelse(V1 == 1, "VIS", 
                           ifelse(V1 == 2, "SMN",
                                  ifelse(V1 == 3, "DAN",
                                         ifelse(V1 == 4, "SAN",
                                                ifelse(V1 == 5, "LIM", 
                                                       ifelse(V1 == 6, "FPN",
                                                               ifelse(V1 == 7, "DMN", "no"))))))))
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

##################### canonical correlation analysis ########################
set.seed(1234)
# filter out health controls
dat_mdd <- dat_td %>% filter(group == "MDD")
plot <- FALSE
grid_start <- 0.0001
grid_end <- 1
grid_length <- 51
dat_cca <- matrix(nrow = 400, ncol = 4, data = NA)
dat_cca_cor <- matrix(nrow = 17, ncol = 4, data = NA)

# regular CCA for saMDD ----------------------------------------------------
dat_amdd <- dat_mdd %>% filter(LPAgroup == "saMDD")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca[,1] <- res.rcc$xcoef[,1]
dat_cca_cor[,1] <- res.rcc$cor
# regular CCA for sdMDD ----------------------------------------------------
dat_amdd <- dat_mdd %>% filter(LPAgroup == "sdMDD")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca[,2] <- res.rcc$xcoef[,1]
dat_cca_cor[,2] <- res.rcc$cor
# regular CCA for miMDD ----------------------------------------------------
dat_amdd <- dat_mdd %>% filter(LPAgroup == "miMDD")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca[,3] <- res.rcc$xcoef[,1]
dat_cca_cor[,3] <- res.rcc$cor
# regular CCA for mdMDD ----------------------------------------------------
dat_amdd <- dat_mdd %>% filter(LPAgroup == "mdMDD")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca[,4] <- res.rcc$xcoef[,1]
dat_cca_cor[,4] <- res.rcc$cor

########################## severe and moderated patients #######################

dat_cca_severe <- matrix(nrow = 400, ncol = 3, data = NA)
dat_cca_severe_cor <- matrix(nrow = 17, ncol = 3, data = NA)

# regular CCA for severe MDD  ------------------------------------------------
dat_amdd <- dat_mdd %>% filter(severity == "severe")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca_severe[,1] <- res.rcc$xcoef[,1]
dat_cca_severe_cor[,1] <- res.rcc$cor[1]

# regular CCA for modulated MDD -----------------------------------------------
dat_amdd <- dat_mdd %>% filter(severity == "moderated")
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca_severe[,2] <- res.rcc$xcoef[,1]
dat_cca_severe_cor[,2] <- res.rcc$cor[1]

# regular CCA for all MDD -----------------------------------------------
dat_amdd <- dat_mdd 
y <- dat_amdd %>% select(starts_with("HAMD_wave1_item"))
x <- dat_amdd %>% select(starts_with("R")) %>% select(-rawID)
res.regular <- CCA::estim.regul(X = x, Y = y, plt = plot,
                                grid1 = seq(grid_start, grid_end, length = grid_length),
                                grid2 = seq(grid_start, grid_end, length = grid_length))

res.rcc <- CCA::rcc(X = x, Y = y, lambda1 = res.regular$lambda1, 
                    lambda2 = res.regular$lambda2)

dat_cca_severe[,3] <- res.rcc$xcoef[,1]
dat_cca_severe_cor[,3] <- res.rcc$cor[1]

# obtain and visualize the CCA loading of the first component
p.cca_severe <- DrawMultipleSurfaceOnSchaefer400(dat_cca_severe, c("severe","modulated","All"),
                                          p_title = "CCA first component", 
                                          legend_title = "loading", highColor = "red")

Cairo::Cairo(file = "outputs/Analysis_2_CCA_severe.png",
             width = 30, height = 20, units = "cm", dpi = 300)
print(p.cca_severe)
dev.off()

############################# collect the CCA results ########################

# for CCA correlation 
colnames(dat_cca_cor) <- c("saMDD","sdMDD","miMDD","mdMDD")
colnames(dat_cca_severe_cor) <- c("severe","modulated","All")

dat_cca_cor_df <- data.frame(t(dat_cca_cor)) %>%
  mutate(group = c("saMDD","sdMDD","miMDD","mdMDD"), .before = 1)
colnames(dat_cca_cor_df) <- paste0("cor",formatC(1:17, width = 2, flag = "0"))

dat_cca_severe_cor_df <- data.frame(t(dat_cca_severe_cor)) %>%
  mutate(group = c("severe","modulated","All"), .before = 1)
colnames(dat_cca_severe_cor_df) <- paste0("cor",formatC(1:17, width = 2, flag = "0"))

dat_cca_cor_all_df <- rbind(dat_cca_cor_df, dat_cca_severe_cor_df)

rio::export(dat_cca_cor_all_df, file = "inputs/Analysis_CCA_correlation.xlsx")

# for CCA loading
colnames(dat_cca) <- c("saMDD","sdMDD","miMDD","mdMDD")
colnames(dat_cca_severe) <- c("severe","modulated","All")

dat_cca_lpa_df <- data.frame(t(dat_cca)) %>%
  mutate(group = c("saMDD","sdMDD","miMDD","mdMDD"), .before = 1)
colnames(dat_cca_lpa_df) <- paste0("R",formatC(1:400, width = 3, flag = "0"))

dat_cca_severe_df <- data.frame(t(dat_cca_severe)) %>%
  mutate(group = c("severe","modulated","All"), .before = 1)
colnames(dat_cca_severe_df) <- paste0("R",formatC(1:400, width = 3, flag = "0"))

dat_cca_df <- rbind(dat_cca_lpa_df, dat_cca_severe_df)

rio::export(dat_cca_df, file = "inputs/Anslysis_CCA_loading.xlsx")

############################# visualization ########################
dat_plt_cca <- rio::import("inputs/Anslysis_CCA_loading.xlsx")
dat_plt_tdp <- rio::import("inputs/Analysis2_TDp_difference.xlsx")

dat_plt_cca_mod <- dat_plt_cca[1:4,2:ncol(dat_plt_cca)]
dat_plt_cca_mod <- t(dat_plt_cca_mod)
dat_plt_cca_mod <- data.frame(dat_plt_cca_mod)
colnames(dat_plt_cca_mod) <- c("saMDD", "sdMDD", "miMDD", "mdMDD")
name_collect <-  c("saMDD", "sdMDD", "miMDD", "mdMDD")

# brain surface plot
p_cca <- DrawMultipleSurfaceOnSchaefer400(dat_plt_cca_mod, name_collect,
                                 p_title = "CCA first component", 
                                 legend_title = "loading", legent_pos = "left")

# CCA between-group correlation
test <- bruceR::Corr(dat_plt_cca_mod, method = "spearman", digits = 3, p.adjust = "fdr")
p_cor_cca <- test$plot

################ coupling between td difference and CCA loading ###############
cor.test(dat_plt_tdp$saMDD, dat_plt_cca_mod$saMDD, method = "spearman") # 0.07472878, p-value = 0.1357
cor.test(dat_plt_tdp$sdMDD, dat_plt_cca_mod$sdMDD, method = "spearman") # -0.03104944, p-value = 0.5356
cor.test(dat_plt_tdp$miMDD, dat_plt_cca_mod$miMDD, method = "spearman") # 0.05582266, p-value = 0.2652
cor.test(dat_plt_tdp$mdMDD, dat_plt_cca_mod$mdMDD, method = "spearman") # -0.03055482, p-value = 0.5421

dat_cca_tdp <- data.frame(subtype = name_collect, corr = c(.074, -.031, .055, -.031))
dat_cca_tdp$subtype <- factor(dat_cca_tdp$subtype, levels = c("saMDD", "sdMDD", "miMDD", "mdMDD"))
p_bar_ccatdp <- ggplot(dat_cca_tdp, aes(x = subtype, y = corr)) +
  geom_bar(stat = "identity", width = .6, fill = "steelblue") +
  geom_text(aes(label = corr), position = position_stack(vjust = 1.1)) +
  geom_hline(yintercept = 0) + ylab("Correlation") +
  theme_classic() +
  easy_text_size(13) + easy_remove_x_axis(what = "title") + easy_rotate_x_labels(angle = -30)

# merge all plots
p_final <- p_cca / (p_cor_cca + p_bar_ccatdp) +
  plot_layout(heights = c(2,1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold"))
p_final
ggsave(plot = p_final, filename = "outputs/Fig3.png", width = 10, height = 10)
ggsave(plot = p_final, filename = "outputs/Fig3.tiff", width = 10, height = 10)
