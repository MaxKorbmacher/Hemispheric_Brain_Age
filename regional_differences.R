# 1) regional differences between hemispheres
# 2) features' age sensitivity
# 30.10.2023 Max Korbmacher (max.korbmacher@gmail.com)

## PREP ####

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

# read in data
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/T1.csv")
dMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/dMRI.csv")
#multimodal = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/combi.csv")
demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")
# HEMISHPERIC DIFFERENCES BOTH SEXES ####
# select the data for each hemisphere and get mean scores for each of the metrics
# T1: surface area, thickness, volume
T1_lh = T1 %>% dplyr::select(starts_with(c("lh", "Left")))
T1_rh = T1 %>% dplyr::select(starts_with(c("rh","Right")))
T1_t.tests = list()
effects = c()
for (i in 1:ncol(T1_lh)){
  x = T1_lh[i]
  y = T1_rh[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  T1_t.tests[[i]] = t.test(df$x, df$y, paired = T)
  effects[i] = cohensD(df$x,df$y, method = "paired")
}

T1_p = c()
for (i in 1:length(T1_t.tests)){
  T1_p[i] = (T1_t.tests[[i]]$p.value)
}
T1res = data.frame(names(T1_lh), T1_p, effects)
# absolute value
sum(ifelse(T1_p > 0.05/length(T1_t.tests), 1, 0))
# relative value
sum(ifelse(T1_p > 0.05/length(T1_t.tests), 1, 0))/length(T1_t.tests)
# effect size assessment
T1res %>% filter(T1_p < 0.05/length(T1_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)))

# dMRI: 28 metrics
dMRI_lh = dMRI %>% dplyr::select(ends_with(c("L")))
dMRI_rh = dMRI %>% dplyr::select(ends_with(c("R")))

dMRI_t.tests = list()
effects = c()
for (i in 1:ncol(dMRI_lh)){
  x = dMRI_lh[i]
  y = dMRI_rh[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  dMRI_t.tests[[i]] = t.test(df$x, df$y, paired = T)
  effects[i] = cohensD(df$x,df$y, method = "paired")
}

dMRI_p = c()
for (i in 1:length(dMRI_t.tests)){
  dMRI_p[i] = (dMRI_t.tests[[i]]$p.value)
}
dMRIres = data.frame(names(dMRI_lh), dMRI_p, effects)
# absolute value
sum(ifelse(dMRI_p > 0.05/length(dMRI_t.tests), 1, 0))
# relative value
sum(ifelse(dMRI_p > 0.05/length(dMRI_t.tests), 1, 0))/length(dMRI_t.tests)
dMRIres %>% filter(dMRI_p < 0.05/length(dMRI_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)))

# SEX STRATIFIED ####
#
# starting with T1-weighted
#
stratdat = merge(demo, T1, by ="eid")
T1_lh_m = stratdat %>% filter(Sex == "Male") %>% dplyr::select(starts_with(c("lh", "Left")))
T1_rh_m = stratdat %>% filter(Sex == "Male") %>% dplyr::select(starts_with(c("rh","Right")))
T1_lh_f = stratdat %>% filter(Sex == "Female") %>% dplyr::select(starts_with(c("lh", "Left")))
T1_rh_f = stratdat %>% filter(Sex == "Female") %>% dplyr::select(starts_with(c("rh","Right")))
l_list = list(T1_lh_m, T1_lh_f)
r_list = list(T1_rh_m, T1_rh_f)
T1res_sex = list()
for (n in 1:2){
  T1_t.tests = list()
  effects = c()
for (i in 1:ncol(T1_lh)){
  x = l_list[[n]][i]
  y = r_list[[n]][i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  T1_t.tests[[i]] = t.test(df$x, df$y, paired = T)
  effects[i] = cohensD(df$x,df$y, method = "paired")
}
  T1_p = c()
for (l in 1:length(T1_t.tests)){
  T1_p[l] = (T1_t.tests[[l]]$p.value)
}
T1res_sex[[n]] = data.frame(names(T1_lh), T1_p, effects)
}
# check results: absolute number of non-sigs
T1res_sex[[1]] %>% filter(T1_p > 0.05/length(T1_t.tests)) %>% nrow()# males
T1res_sex[[2]] %>% filter(T1_p > 0.05/length(T1_t.tests)) %>% nrow()# females
# relative value
T1res_sex[[1]] %>% filter(T1_p > 0.05/length(T1_t.tests)) %>% nrow()/length(T1_t.tests)# males
T1res_sex[[2]] %>% filter(T1_p > 0.05/length(T1_t.tests)) %>% nrow()/length(T1_t.tests)# females
# effect size assessment
T1res_sex[[1]] %>% filter(T1_p < 0.05/length(T1_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)), MAX =max(effects), MIN = min(effects))# males
T1res_sex[[2]] %>% filter(T1_p < 0.05/length(T1_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)), MAX =max(effects), MIN = min(effects))# females
# make table
sex_T1_diff = merge(T1res_sex[1], T1res_sex[2], by = "names.T1_lh.")
names(sex_T1_diff) = c("Region", "p_male","Cohensd_male", "p_female","Cohensd_female")
write.csv(sex_T1_diff, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_sex_T1_diff.csv")
#
#
# diffusion data
stratdat = merge(demo, dMRI, by ="eid")
dMRI_lh_m = stratdat %>% filter(Sex == "Male") %>% dplyr::select(names(dMRI_lh))
dMRI_rh_m = stratdat %>% filter(Sex == "Male") %>% dplyr::select(names(dMRI_rh))
dMRI_lh_f = stratdat %>% filter(Sex == "Female") %>% dplyr::select(names(dMRI_lh))
dMRI_rh_f = stratdat %>% filter(Sex == "Female") %>% dplyr::select(names(dMRI_rh))
l_list = list(dMRI_lh_m, dMRI_lh_f)
r_list = list(dMRI_rh_m, dMRI_rh_f)
dMRIres_sex = list()
for (n in 1:2){
  dMRI_t.tests = list()
  effects = c()
  for (i in 1:ncol(dMRI_lh)){
    x = l_list[[n]][i]
    y = r_list[[n]][i]
    df = data.frame(x,y)
    colnames(df) = c("x","y")
    dMRI_t.tests[[i]] = t.test(df$x, df$y, paired = T)
    effects[i] = cohensD(df$x,df$y, method = "paired")
  }
  dMRI_p = c()
  for (l in 1:length(dMRI_t.tests)){
    dMRI_p[l] = (dMRI_t.tests[[l]]$p.value)
  }
  dMRIres_sex[[n]] = data.frame(names(dMRI_lh), dMRI_p, effects)
}
# check results: absolute number of non-sigs
dMRIres_sex[[1]] %>% filter(dMRI_p > 0.05/length(dMRI_t.tests)) %>% nrow()# males
dMRIres_sex[[2]] %>% filter(dMRI_p > 0.05/length(dMRI_t.tests)) %>% nrow()# females
# relative value
dMRIres_sex[[1]] %>% filter(dMRI_p > 0.05/length(dMRI_t.tests)) %>% nrow()/length(dMRI_t.tests)# males
dMRIres_sex[[2]] %>% filter(dMRI_p > 0.05/length(dMRI_t.tests)) %>% nrow()/length(dMRI_t.tests)# females
# effect size assessment
dMRIres_sex[[1]] %>% filter(dMRI_p < 0.05/length(dMRI_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)), MAX =max(effects), MIN = min(effects))# males
dMRIres_sex[[2]] %>% filter(dMRI_p < 0.05/length(dMRI_t.tests)) %>% summarize(ES = mean(abs(effects)), SD = sd(abs(effects)), MAX =max(effects), MIN = min(effects))# females
# make table
sex_dMRI_diff = merge(dMRIres_sex[1], dMRIres_sex[2], by = "names.dMRI_lh.")
names(sex_dMRI_diff) = c("Region", "p_male","Cohensd_male", "p_female","Cohensd_female")
write.csv(sex_dMRI_diff, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_sex_dMRI_diff.csv")

# AGE-SENSITIVITY BOTH SEXES ####
## T1-weighted ####
stratdat = merge(demo, T1, by ="eid")
T1_LRT = stratdat %>% dplyr::select(Age, Sex, Scanner, starts_with(c("lh", "Left","rh","Right")))
metrics_list = names(T1_LRT[4:ncol(T1_LRT)])
mod0list = list()
mod1list = list()
SS = c()
F.stat = c()
p = c()

### linear modelling #### 
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~",i, "+Sex+Scanner"))
  mod0list[[i]] = lm(f0,data = T1_LRT)
  mod1list[[i]] = lm(f1,data = T1_LRT)
  l = anova(mod0list[[i]], mod1list[[i]])
  SS[i] = l$`Sum of Sq`[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p_adj = length(p)*p
# make df and save
linear_hemi_effects = data.frame(metrics_list, SS, F.stat, p, p_adj)
write.csv(linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_linear_hemi_effects.csv")
# round to 2 digits for simpler reading
linear_hemi_effects[2:4] = round(linear_hemi_effects[2:4], digits = 2)
# absolute value
paste(sum(ifelse(p_adj > 0.05/length(p_adj), 1, 0)),"of", length(p_adj), "were non-significantly age-related using linear associations")
# relative value
paste("Hence, proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/length(p), 1, 0))/length(p), "percent.")
# effect size assessment
linear_hemi_effects %>% filter(p_adj < 0.05/length(p)) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat)))

### non-linear modelling ####
Deviance = c()
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~s(",i, ",k=4)+Sex+Scanner"))
  mod0list[[i]] = gam(f0,data = T1_LRT, method = "REML")
  mod1list[[i]] = gam(f1,data = T1_LRT, method = "REML")
  l = anova.gam(mod0list[[i]], mod1list[[i]], test = "F")
  Deviance[i] = l$Deviance[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p_adj = length(p)*p
# make df and save
non_linear_hemi_effects = data.frame(metrics_list, Deviance, F.stat, p, p_adj)
write.csv(non_linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_T1_non_linear_hemi_effects.csv")
non_linear_hemi_effects[2:4] = round(non_linear_hemi_effects[2:4], digits = 2)
# absolute value
paste(sum(ifelse(p_adj > 0.05/length(p_adj), 1, 0)),"of", length(p_adj), "were non-significantly age-related using non-linear associations")
# relative value
paste("Hence, proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/length(p), 1, 0))/length(p), "percent.")
# effect size assessment
non_linear_hemi_effects %>% filter(p_adj < 0.05/length(p)) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat)))



## Diffusion-weighted ####
stratdat = merge(demo, dMRI, by ="eid")
dMRI_LRT = stratdat %>% dplyr::select(Age, Sex, Scanner, ends_with(c("L","R")))
metrics_list = names(dMRI_LRT[4:ncol(dMRI_LRT)])
mod0list = list()
mod1list = list()
SS = c()
F.stat = c()
p = c()

### linear modelling #### 
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~",i, "+Sex+Scanner"))
  mod0list[[i]] = lm(f0,data = dMRI_LRT)
  mod1list[[i]] = lm(f1,data = dMRI_LRT)
  l = anova(mod0list[[i]], mod1list[[i]])
  SS[i] = l$`Sum of Sq`[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p_adj = length(p)*p
# make df and save
linear_hemi_effects = data.frame(metrics_list, SS, F.stat, p, p_adj)
write.csv(linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_linear_hemi_effects.csv")
# round to 2 digits for simpler reading
linear_hemi_effects[2:4] = round(linear_hemi_effects[2:4], digits = 2)
# absolute value
paste(sum(ifelse(p_adj > 0.05/length(p_adj), 1, 0)),"of", length(p_adj), "were non-significantly age-related using linear associations")
# relative value
paste("Hence, proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/length(p), 1, 0))/length(p), "percent.")
# effect size assessment
linear_hemi_effects %>% filter(p_adj < 0.05/length(p)) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat)))

### non-linear modelling ####
Deviance = c()
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~s(",i, ",k=4)+Sex+Scanner"))
  mod0list = gam(f0,data = dMRI_LRT, method = "REML")
  mod1list = gam(f1,data = dMRI_LRT, method = "REML")
  l = anova.gam(mod0list, mod1list, test = "F")
  Deviance[i] = l$Deviance[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p_adj = length(p)*p
# make df and save
non_linear_hemi_effects = data.frame(metrics_list, Deviance, F.stat, p, p_adj)
write.csv(non_linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_dMRI_non_linear_hemi_effects.csv")
non_linear_hemi_effects[2:4] = round(non_linear_hemi_effects[2:4], digits = 2)
# absolute value
paste(sum(ifelse(p_adj > 0.05/length(p_adj), 1, 0)),"of", length(p_adj), "were non-significantly age-related using non-linear associations")
# relative value
paste("Hence, proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/length(p), 1, 0))/length(p), "percent.")
# effect size assessment
non_linear_hemi_effects %>% filter(p_adj < 0.05/length(p)) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat)))