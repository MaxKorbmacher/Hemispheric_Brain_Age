# Sex stratified regional age-sensitivity
# due to high computational costs, this was written out as stand-alone script
# 01.11.2023, Max Korbmacher (max.korbmacher@gmail.com)
#
# PREP ####
# packages
library(tidyquant)
library(afex)
library(effects)
library(dplyr)
library(reshape)
library(ggplot2)
library(tidyr)
library(viridis)
library(lme4)
library(MuMIn)
library(effectsize)
library(ggpubr)
library(cowplot)
library(ggplotify)
library(pheatmap)
library(RColorBrewer)
library(ppcor)
library(MuMIn)
library(ggdist)
library(emmeans)
library(dataPreparation)
library(lsr)
library(mgcv)
#
# read in data
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/T1.csv")
dMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/dMRI.csv")
#multimodal = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/combi.csv")
demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")
#
# COMPUTE ####
# define empty components to be filled
# mod0list = list()
# mod1list = list()
#
# DIFFUSION ####
stratdat = merge(demo, dMRI, by ="eid")
dMRI_LRT_m = stratdat %>% dplyr::filter(Sex == "Male") %>% dplyr::select(Age, Scanner, ends_with(c("L","R")))
dMRI_LRT_f = stratdat %>% dplyr::filter(Sex == "Female") %>% dplyr::select(Age, Scanner, ends_with(c("L","R")))
metrics_list = names(dMRI_LRT_m[3:ncol(dMRI_LRT_m)])
dMRI_LRT = list(dMRI_LRT_m, dMRI_LRT_f)
#
### linear modeling #### 
print("Starting linear modelling.")
linear_hemi_effects = list()
for (df in 1:length(dMRI_LRT)){
  dat = dMRI_LRT[[df]]
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  for (i in (metrics_list)){
  f0 = formula(paste("Age~Scanner"))
  f1 = formula(paste("Age~",i, "+Scanner"))
  mod0list = lm(f0,data = dat)
  mod1list = lm(f1,data = dat)
  l = anova(mod0list, mod1list)
  SS[i] = l$`Sum of Sq`[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
  print(paste("LRT assessing the LINEAR age-sensitivity of", i, "completed."))
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list, SS, F.stat, p, p_adj)
  print(paste("Data set ", df, " processed."))
}
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_linear_hemi_effects_MALES.csv")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_linear_hemi_effects_FEMALES.csv")
#
### non-linear modeling ####
non_linear_hemi_effects = list()
print("Starting non-linear modelling.")
for (df in 1:length(dMRI_LRT)){
  dat = dMRI_LRT[[df]]
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  for (i in (metrics_list)){
  f0 = formula(paste("Age~Scanner"))
  f1 = formula(paste("Age~s(",i, ",k=4)+Scanner"))
  mod0list = gam(f0,data = dat, method = "REML")
  mod1list = gam(f1,data = dat, method = "REML")
  l = anova.gam(mod0list, mod1list, test = "F")
  Deviance[i] = l$Deviance[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
  paste("LRT assessing the age-sensitivity of NON-LINEAR", i, "completed.")
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  non_linear_hemi_effects[[df]] = data.frame(metrics_list, Deviance, F.stat, p, p_adj)
  paste("Data set ", df, " processed.")
}
write.csv(non_linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_non_linear_hemi_effects_MALES.csv")
write.csv(non_linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_non_linear_hemi_effects_FEMALES.csv")
# T1 ####
stratdat = merge(demo, T1, by ="eid")
T1_LRT_m = stratdat %>% dplyr::filter(Sex == "Male") %>% dplyr::select(Age, Scanner, starts_with(c("lh", "Left","rh","Right")))
T1_LRT_f = stratdat %>% dplyr::filter(Sex == "Female") %>% dplyr::select(Age, Scanner, starts_with(c("lh", "Left","rh","Right")))
metrics_list = names(T1_LRT_m[3:ncol(T1_LRT_m)])
T1_LRT = list(T1_LRT_m, T1_LRT_f)
#
### linear modeling #### 
linear_hemi_effects = list()
for (df in 1:length(T1_LRT)){
  dat = T1_LRT[[df]]
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  for (i in (metrics_list)){
    f0 = formula(paste("Age~Scanner"))
    f1 = formula(paste("Age~",i, "+Scanner"))
    mod0list = lm(f0,data = dat)
    mod1list = lm(f1,data = dat)
    l = anova(mod0list, mod1list)
    SS[i] = l$`Sum of Sq`[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    paste("LRT assessing the age-sensitivity of", i, "completed.")
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list, SS, F.stat, p, p_adj)
  paste("Data set ", df, " processed.")
}
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_linear_hemi_effects_MALES.csv")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_linear_hemi_effects_FEMALES.csv")
#
### non-linear modeling ####
non_linear_hemi_effects = list()
for (df in 1:length(T1_LRT)){
  dat = T1_LRT[[df]]
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  for (i in (metrics_list)){
    f0 = formula(paste("Age~Scanner"))
    f1 = formula(paste("Age~s(",i, ",k=4)+Scanner"))
    mod0list = gam(f0,data = dat, method = "REML")
    mod1list = gam(f1,data = dat, method = "REML")
    l = anova.gam(mod0list, mod1list, test = "F")
    Deviance[i] = l$Deviance[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    print(paste("LRT assessing the age-sensitivity of", i, "completed."))
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  non_linear_hemi_effects[[df]] = data.frame(metrics_list, Deviance, F.stat, p, p_adj)
  paste("Data set ", df, " processed.")
}
write.csv(non_linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_non_linear_hemi_effects_MALES.csv")
write.csv(non_linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_non_linear_hemi_effects_FEMALES.csv")
#
print("Done. Now check the resulting data frames (Hemi_*FEMALES.csv and Hemi_*MALES.csv) in the export folder.")
print("...")
print("some more meta stats on the tables")
#
#
#
#
# reported stats in the paper >> change loaded data.
# dMRI
d1 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_linear_hemi_effects_MALES.csv")
d2 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_linear_hemi_effects_FEMALES.csv")
d3 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_non_linear_hemi_effects_MALES.csv")
d4 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_non_linear_hemi_effects_FEMALES.csv")
# T1
t1 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_linear_hemi_effects_MALES.csv")
t2 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_linear_hemi_effects_FEMALES.csv")
t3 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_non_linear_hemi_effects_MALES.csv")
t4 = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_non_linear_hemi_effects_FEMALES.csv")
#non_linear_hemi_effects[2:4] = round(non_linear_hemi_effects[2:4], digits = 2)
res = list(d1, d2, d3, d4, t1, t2, t3, t4)
ress = c("dMRI lin males", "dMRI lin females", "dMRI NON-lin males", "dMRI NON-lin females", "t1 lin males", "t1 lin females", "t1 NON-lin males", "t1 NON-lin females")
# absolute value
for (i in 1:length(res)){
  print(ress[i])
  p_adj = res[[i]][ncol(res[[i]])]
  print("Total number of non-sig/non-age sens features:")
  print(paste(nrow(res[[i]])-sum(ifelse(p_adj > 0.05/nrow(res[[i]]), 1, 0)),"of", nrow(res[[i]]), "were non-significantly age-related using non-linear associations"))
  print("Percentages:")
  print(paste("Proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/nrow(res[[i]]), 1, 0))/nrow(res[[i]]), "percent."))
  print("mean absolute F-statistic (SD)")
  print(res[[i]] %>% filter(p_adj < 0.05/nrow(res[[i]])) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat))))
}
