##########################
#
# Subjects information description
#
# This script is used to summary 
# the demographic information about 
# subjects between MDD and HC
#
# Liang Qunjun 2023/07/19

library(tidyverse)
library(NbClust)
library(ggiraphExtra)
library(ggsci)
library(ggeasy)
library(tidyLPA)
library(psych)
library(mclust)
library(ggsignif)
library(patchwork)
library(scales)
source("scripts/R/function_PvalueForTable1.R")
source("scripts/R/function_ColorFroGroup.R")

# load the data
dat_sbj <- rio::import('inputs/subject_merge_MDD_HC.xlsx')

# fix the education 
dat_sbj <- dat_sbj %>% mutate(educations = ifelse(education == 1, 'Illiterate',
                                                  ifelse(education == 2, 'Primary education',
                                                         ifelse(education == 3, 'Junior high school',
                                                                ifelse(education == 4, 'Senior high school',
                                                                       ifelse(education == 5, 'Undergraduate', 'Graduate'))))))

# summary the clinical scales
dat_mdd <- dat_sbj %>% filter(group == "MDD") 
dat_hc <- dat_sbj %>% filter(group == "HC")

#############################################################
#
# calculate the factor and total score
#
#############################################################

# define the items
item_depression <- c(1,2,7,8,13)
item_anxiety <- c(15,9,10,11,16)
item_sleepness <- c(4,5,6)

item_A <- c(1,2,7,8,10,13)
item_B <- c(4,5,6,9,11,12,14,15,17)
item_C <- c(3,16)

# calculate the factor and total score
hamd_wave1 <- dat_mdd %>% select(participant_id, contains('HAMD_wave1_item')) 
hamd_wave1 <- hamd_wave1 %>%
  mutate(HAMD_wave1_total = rowSums(hamd_wave1[,2:ncol(hamd_wave1)])) %>%
  mutate(HAMD_wave1_depression = rowSums(hamd_wave1[, 1+item_depression])) %>%
  mutate(HAMD_wave1_anxiety = rowSums(hamd_wave1[, 1+item_anxiety])) %>%
  mutate(HAMD_wave1_sleepness = rowSums(hamd_wave1[, 1+item_sleepness])) %>%
  mutate(HAMD_wave1_A = rowSums(hamd_wave1[, 1+item_A])) %>%
  mutate(HAMD_wave1_B = rowSums(hamd_wave1[, 1+item_B])) %>%
  mutate(HAMD_wave1_C = rowSums(hamd_wave1[, 1+item_C]))

hamd_wave2 <- dat_mdd %>% select(participant_id, contains('HAMD_wave2_item')) 
hamd_wave2 <- hamd_wave2 %>%
  mutate(HAMD_wave2_total = rowSums(hamd_wave2[,2:ncol(hamd_wave2)])) %>%
  mutate(HAMD_wave2_depression = rowSums(hamd_wave2[, 1+item_depression])) %>%
  mutate(HAMD_wave2_anxiety = rowSums(hamd_wave2[, 1+item_anxiety])) %>%
  mutate(HAMD_wave2_sleepness = rowSums(hamd_wave2[, 1+item_sleepness])) %>%
  mutate(HAMD_wave2_A = rowSums(hamd_wave2[, 1+item_A])) %>%
  mutate(HAMD_wave2_B = rowSums(hamd_wave2[, 1+item_B])) %>%
  mutate(HAMD_wave2_C = rowSums(hamd_wave2[, 1+item_C]))

hamd_wave3 <- dat_mdd %>% select(participant_id, contains('HAMD_wave3_item')) 
hamd_wave3 <- hamd_wave3 %>%
  mutate(HAMD_wave3_total = rowSums(hamd_wave3[,2:ncol(hamd_wave3)])) %>%
  mutate(HAMD_wave3_depression = rowSums(hamd_wave3[, 1+item_depression])) %>%
  mutate(HAMD_wave3_anxiety = rowSums(hamd_wave3[, 1+item_anxiety])) %>%
  mutate(HAMD_wave3_sleepness = rowSums(hamd_wave3[, 1+item_sleepness])) %>%
  mutate(HAMD_wave3_A = rowSums(hamd_wave3[, 1+item_A])) %>%
  mutate(HAMD_wave3_B = rowSums(hamd_wave3[, 1+item_B])) %>%
  mutate(HAMD_wave3_C = rowSums(hamd_wave3[, 1+item_C]))

# calculate HAMA wave1 total score
hama_wave1 <- dat_mdd %>% select(participant_id, contains('HAMA_wave1_item')) 
hama_wave1 <- hama_wave1 %>% mutate(HAMA_wave1_total = rowSums(hama_wave1[,2:ncol(hama_wave1)]) )

# merger the factor to the main data table
dat_mdd_factor <- dat_mdd %>% left_join(hamd_wave1) %>% left_join(hamd_wave2) %>%
  left_join(hamd_wave3) %>% left_join(hama_wave1)

dat_mdd_factor <- dat_mdd_factor %>% filter(HAMD_wave1_total >= 17)

##############################################################
#
# Latent profile analysis - LPA
#
##############################################################
set.seed(1234)
dat_mdd_factor %>% select(starts_with("HAMD_wave1_item")) %>%
  select(1:17) %>% 
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1) %>%
  compare_solutions(statistics = c('AIC','BIC',"AWE", "CLC", "KIC"))
# Model Classes AIC      BIC      AWE      CLC      KIC     
# 1     2       5368.773 5519.847 5928.938 5266.756 5423.773
# 1     3       5414.744 5618.113 6169.619 5276.608 5487.744
# 1     4       5187.368 5443.032 6136.784 5013.280 5278.368
# 
# Best model according to AIC is Model 1 with 4 classes.
# Best model according to BIC is Model 1 with 4 classes.
# Best model according to AWE is Model 1 with 2 classes.
# Best model according to CLC is Model 1 with 4 classes.
# Best model according to KIC is Model 1 with 4 classes.
# 
# An analytic hierarchy process, 
# based on the fit indices AIC, AWE, BIC, CLC, and KIC (Akogul & Erisoglu, 2017), 
# suggests the best solution is Model 1 with 4 classes.

# call mixture model
model_lpa <- dat_mdd_factor %>% select(starts_with("HAMD_wave1_item")) %>%
  estimate_profiles(n_profiles = 2:4, eplot_profiles = 1)

# plot the result
p_lpa_line <- model_lpa$model_1_class_4$estimates %>%
  filter(Category == "Means") %>% 
  mutate(item = rep(paste0("item", formatC(1:17,width = 2, flag = "0")), times = 4)) %>%
  mutate(LPAgroup = ifelse(Class == 1, "saMDD", 
                           ifelse(Class == 2, "miMDD",
                                  ifelse(Class == 3, "mdMDD", "sdMDD")))) %>%
  ggplot(aes(x = item, y = Estimate, group = LPAgroup, color = LPAgroup)) +
  geom_line() + 
  geom_point() +
  scale_color_manual(values=c(ColorForGroup("mdMDD"),ColorForGroup("miMDD"), ColorForGroup("saMDD"), ColorForGroup("sdMDD"))) +
  xlab("HAMD-17 items") +
  theme_classic() + easy_remove_x_axis(what = "title") +
  easy_text_size(15) + easy_rotate_x_labels(angle = -45)
p_lpa_line

Cairo::Cairo(file = 'outputs/Anay1_LPA_lineplot.png')
print(p_lpa_line)
dev.off()
##########################################################
#
# Merge the LPA group and testing HAMD and HAMA difference
#
##########################################################
## merge the data
dat_mdd_factor_lpa <- model_lpa$model_1_class_4$dff %>% 
  mutate(LPAgroup = ifelse(Class == 1, "saMDD", 
                           ifelse(Class == 2, "miMDD",
                                  ifelse(Class == 3, "mdMDD", "sdMDD")))) %>%
  mutate(LPAgroup = factor(LPAgroup, c("saMDD","sdMDD","mdMDD","miMDD"))) %>%
  mutate(severity = ifelse(LPAgroup == "mdMDD" | LPAgroup == "miMDD", "moderated", "severe")) %>%
  select(Class, LPAgroup, severity) %>% 
  cbind(dat_mdd_factor)

## Testing difference in HAMD -----------------------------------------
bruceR::MANOVA(
  data = dat_mdd_factor_lpa,
  subID = "participant_id",
  between = "LPAgroup",
  dv = "HAMD_wave1_total",
  covariate = c("age","gender")
) %>%
  bruceR::EMMEANS(effect = "LPAgroup", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df       t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────
# sdMDD - saMDD   -1.794 (0.936) 129  -1.918  .069 .   -0.471 [-1.128,  0.187]
# mdMDD - saMDD   -9.702 (1.010) 129  -9.602 <.001 *** -2.545 [-3.255, -1.835]
# mdMDD - sdMDD   -7.907 (0.870) 129  -9.088 <.001 *** -2.074 [-2.686, -1.463]
# miMDD - saMDD  -11.458 (1.031) 129 -11.110 <.001 *** -3.006 [-3.731, -2.281]
# miMDD - sdMDD   -9.664 (0.933) 129 -10.355 <.001 *** -2.535 [-3.191, -1.879]
# miMDD - mdMDD   -1.756 (1.009) 129  -1.742  .084 .   -0.461 [-1.170,  0.248]
# ─────────────────────────────────────────────────────────────────────────────

## boxlplot to show the HAMD wave1 difference
p_lpa_hamd_diff <- ggplot(dat_mdd_factor_lpa, aes(x = LPAgroup, y = HAMD_wave1_total, fill = LPAgroup)) +
  geom_violin(alpha = .4, color = "white") +
  geom_boxplot(alpha = .5, width = .3) +
  geom_signif(annotations = rep("***", times = 4), 
              textsize = 7, vjust = .7,
              y_position = c(46.5, 45, 43.5, 42), 
              xmin = c(2, 1, 2, 1), 
              xmax = c(4, 4, 3, 3),
              tip_length = 0) +
  scale_fill_lancet() + ylab("HAMD score") + xlab("LPA subgroups") +
  scale_fill_manual(values=c(ColorForGroup("saMDD"),ColorForGroup("sdMDD"), ColorForGroup("mdMDD"), ColorForGroup("miMDD"))) +
  theme_classic() +
  easy_text_size(15) + easy_remove_x_axis(what = "title")
p_lpa_hamd_diff

## HAMA difference --------------------------------
## MANOVA with bruceR
bruceR::MANOVA(
  data = dat_mdd_factor_lpa,
  subID = "participant_id",
  between = "LPAgroup",
  dv = "HAMA_wave1_total",
  covariate = c("age","gender")
) %>%
  bruceR::EMMEANS(effect = "LPAgroup", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────
# sdMDD - saMDD   -3.078 (1.463) 127 -2.105  .045 *   -0.522 [-1.186,  0.143]
# mdMDD - saMDD   -7.066 (1.593) 127 -4.437 <.001 *** -1.197 [-1.921, -0.474]
# mdMDD - sdMDD   -3.987 (1.360) 127 -2.932  .006 **  -0.676 [-1.293, -0.058]
# miMDD - saMDD   -8.062 (1.612) 127 -5.001 <.001 *** -1.366 [-2.098, -0.634]
# miMDD - sdMDD   -4.984 (1.446) 127 -3.448  .002 **  -0.844 [-1.501, -0.188]
# miMDD - mdMDD   -0.996 (1.578) 127 -0.632  .529     -0.169 [-0.885,  0.548]
# ────────────────────────────────────────────────────────────────────────────

p_lpa_hama_diff <- ggplot(dat_mdd_factor_lpa, aes(x = LPAgroup, y = HAMA_wave1_total, fill = LPAgroup)) +
  geom_violin(alpha = .4, color = "white") +
  geom_boxplot(alpha = .5, width = .3) +
  geom_signif(annotations = c("**","**","***","***","*"), 
              textsize = 7, vjust = .7,
              y_position = c(39, 41, 43, 45, 47), 
              xmin = c(1, 1, 1, 2, 2), 
              xmax = c(2, 3, 4, 3, 4),
              tip_length = 0) +
  scale_fill_manual(values=c(ColorForGroup("saMDD"),ColorForGroup("sdMDD"), ColorForGroup("mdMDD"), ColorForGroup("miMDD"))) +
  ylab("HAMA score") + xlab("LPA subgroups") +
  theme_classic() +
  easy_text_size(15) + easy_remove_x_axis(what = "title")
p_lpa_hama_diff


##############################################################
#
# Merge the plots
#
##############################################################

p_final <- p_lpa_line / (p_lpa_hamd_diff + p_lpa_hama_diff) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(axis.text = element_text(size = 13), plot.tag = element_text(size = 13, face = "bold"))
p_final

ggsave(plot = p_final, filename = "outputs/Fig1.png", width = 8.5, height = 7.5)
ggsave(plot = p_final, filename = "outputs/Fig1.tiff", width = 8.5, height = 7.5)

##############################################################
#
# Export the new subject information
#
##############################################################

### summary the LPA group 
table1::table1(data = dat_mdd_factor_lpa, 
               ~ age + gender + HAMD_wave1_total + 
                 HAMA_wave1_total | LPAgroup)
# merge mdd and HC
dat_merge <- data.table::rbindlist(list(dat_mdd_factor_lpa, dat_hc), fill = T)
dat_merge$LPAgroup <- as.character(dat_merge$LPAgroup)
dat_merge <- dat_merge %>% 
  mutate(subtype_lpa=ifelse(group == 'MDD', LPAgroup, group)) %>%
  mutate(subtype_sever=ifelse(group == "MDD", severity, group))
table1::table1(data = dat_merge, 
               ~ age + gender| group)
rio::export(dat_merge, file = "inputs/Analysis1_subject_table.xlsx")
