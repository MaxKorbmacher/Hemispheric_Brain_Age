#!/usr/bin/env Rscript
# Linear and non-linear modelling of hemispheric average values
# This code corresponds to the hemispheric differences and age sensitivity sections (with a focus on hemispheric averages)
# 30.10.2023, Max Korbmacher (max.korbmacher@gmail.com)
#
#
# PREP ####
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize, gridExtra,ggsignif,rstatix,data.table,ggridges,dplyr, rms, mgcv, marginaleffects, dataPreparation)
#
# set work dir and load data
setwd("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/analyses/")
# 9 SD over mean included data:
orig_df = read.csv("orig_df.csv")
# save raw data frame (before scaling)
raw = orig_df
# 9 SD over mean excluded data:
# plot_df = remove_sd_outlier(orig_df, n_sigmas = 9) # or load:
plot_df = read.csv("hemispheric_plot_df.csv")
#
# make a list of all the metric names
metrics_list = names(plot_df)[6:ncol(plot_df)]#
#
# SCALE data frames
plot_df[,6:ncol(plot_df)] = plot_df[,6:ncol(plot_df)]%>%mutate_if(is.numeric,scale)
orig_df[,6:ncol(orig_df)] = orig_df[,6:ncol(orig_df)]%>%mutate_if(is.numeric,scale)
#
# relocate some values
plot_df = plot_df %>% dplyr::relocate(X, eid, Age, Sex, Scanner, .after = last_col())
orig_df = orig_df %>% dplyr::relocate(X, eid, Age, Sex, Scanner, .after = last_col())
#
#
# and a list of print version names
betternames = c("BRIA - V intra",      "BRIA - V extra"  ,    "BRIA - V CSF"     ,   "BRIA - microRD"   ,  "BRIA - microFA"   ,  "BRIA - microAX"  ,   "BRIA - microADC"  , 
                "BRIA - dradextra" ,  "BRIA - daxintra"  ,  "BRIA - daxextra"  ,  "DKI - AK"     ,      "DKI - RK"     ,      "DKI - MK"     ,      "DTI - AD"  ,        
                "DTI - RD"    ,       "DTI - MD"    ,       "DTI - FA"   ,        "SMT - FA"     ,      "SMT - MD"       ,    "SMT - trans"    ,    "SMT - long"   ,     
                "SMTmc - diff"      ,    "SMTmc - extraMD" ,   "SMTmc - extra trans", "SMTmc - intra"  ,    "WMTI - AWF"     ,    "WMTI - radEAD"    ,  "WMTI - axEAD"      ,
                "T1 - thickness"   ,  "T1 - area"  ,        "T1 - volume")
betternames_per_hemi = c("BRIA - V intra left","BRIA - V intra right",
                         "BRIA - V extra left"  ,"BRIA - V extra right"  ,
                         "BRIA - V CSF left"     ,"BRIA - V CSF right"     ,
                         "BRIA - microRD left"   ,  "BRIA - microRD right"   ,
                         "BRIA - microFA left"   , "BRIA - microFA right"   ,
                         "BRIA - microAX left"  ,"BRIA - microAX right"  ,
                         "BRIA - microADC left"  , "BRIA - microADC right"  ,
                         "BRIA - dradextra left" ,  "BRIA - dradextra right" , 
                         "BRIA - daxintra left"  ,  "BRIA - daxintra right"  , 
                         "BRIA - daxextra left"  ,  "BRIA - daxextra right"  ,  
                         "DKI - AK left"     ,"DKI - AK right"     ,
                         "DKI - RK left"     ,"DKI - RK right"     ,
                         "DKI - MK left"     ,"DKI - MK right"     ,
                         "DTI - AD left"  ,"DTI - AD right"  ,
                         "DTI - RD left"    ,"DTI - RD right"    ,
                         "DTI - MD left"    ,"DTI - MD right"    ,
                         "DTI - FA left"   ,"DTI - FA right"   ,
                         "SMT - FA left"     ,"SMT - FA right"     ,
                         "SMT - MD left"       ,"SMT - MD right"       ,
                         "SMT - trans left"    ,"SMT - trans right"    ,
                         "SMT - long left"   ,    "SMT - long right"   ,  
                         "SMTmc - diff left","SMTmc - diff right",
                         "SMTmc - extraMD left" ,"SMTmc - extraMD right" ,
                         "SMTmc - extra trans left", "SMTmc - extra trans right", 
                         "SMTmc - intra left"  ,"SMTmc - intra right"  ,
                         "WMTI - AWF left"     ,"WMTI - AWF right"     ,
                         "WMTI - radEAD left"    ,"WMTI - radEAD right"    ,
                         "WMTI - axEAD left"      ,"WMTI - axEAD right"      ,
                         "T1 - thickness left"   ,"T1 - thickness right"   ,
                         "T1 - area left"  ,"T1 - area right"  ,
                         "T1 - volume left", "T1 - volume right")
#
# also create a new data frame >> long format for hemisphere (with dummy indicating hemisphere: L or R)
LEFT = orig_df %>% dplyr::select(contains("lh"), Age, Sex, Scanner)
RIGHT = orig_df %>% dplyr::select(contains("rh"), Age, Sex, Scanner)
some_metrics = c("BRIA_vintra",      "BRIA_vextra"  ,    "BRIA_vcsf"     ,   "BRIA_micrord"   ,  "BRIA_microfa"   ,  "BRIA_microax"  ,   "BRIA_microadc"  , 
                 "BRIA_dradextra" ,  "BRIA_daxintra"  ,  "BRIA_daxextra"  ,  "DKI_AK"     ,      "DKI_RK"     ,      "DKI_MK"     ,      "DTI_AD"  ,        
                 "DTI_RD"    ,       "DTI_MD"    ,       "DTI_FA"   ,        "SMT_FA"     ,      "SMT_MD"       ,    "SMT_trans"    ,    "SMT_long"   ,     
                 "SMTmc_d"      ,    "SMTmc_extramd" ,   "SMTmc_extratrans", "SMTmc_intra"  ,    "WMTI_awf"     ,    "WMTI_radead"    ,  "WMTI_axead"      ,
                 "T1_thickness"   ,  "T1_area"  ,        "T1_volume")
names(LEFT) = c("BRIA_vintra",      "BRIA_vextra"  ,    "BRIA_vcsf"     ,   "BRIA_micrord"   ,  "BRIA_microfa"   ,  "BRIA_microax"  ,   "BRIA_microadc"  , 
                "BRIA_dradextra" ,  "BRIA_daxintra"  ,  "BRIA_daxextra"  ,  "DKI_AK"     ,      "DKI_RK"     ,      "DKI_MK"     ,      "DTI_AD"  ,        
                "DTI_RD"    ,       "DTI_MD"    ,       "DTI_FA"   ,        "SMT_FA"     ,      "SMT_MD"       ,    "SMT_trans"    ,    "SMT_long"   ,     
                "SMTmc_d"      ,    "SMTmc_extramd" ,   "SMTmc_extratrans", "SMTmc_intra"  ,    "WMTI_awf"     ,    "WMTI_radead"    ,  "WMTI_axead"      ,
                "T1_thickness"   ,  "T1_area"  ,        "T1_volume"     ,   "Age", "Sex", "Scanner")
names(RIGHT) = names(LEFT)
both = rbind(LEFT, RIGHT)
both$Hemisphere = c(replicate(nrow(LEFT),"L"),replicate(nrow(RIGHT),"R"))
both$Hemisphere = factor(both$Hemisphere)
write.csv(both, "both.csv")
#
################# LINEAR PLOTTING ####
# ADJUSTED ######
# estimate lms for global scores and plot them
plot_list = list()
for (i in 1:(length(metrics_list))){
  # starting with linear fit
  formula1 = paste(metrics_list[i], "~ Age+Sex+Scanner")
  yname=betternames_per_hemi[i]
  mod = lm(formula1, orig_df)
  plot_list[[i]] = plot_predictions(mod, condition = "Age", rug = FALSE, draw = TRUE, vcov = TRUE) + 
    theme_bw()+ theme(legend.position = "none") + labs(y = yname)
}
# plot for tsd
age_curves = ggarrange(plotlist = plot_list) %>% 
  ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_NEW_LINEAR.pdf", plot = ., width = 15, height = 25, device = "pdf")
# plot for export
age_curves = ggarrange(plotlist = plot_list) %>% 
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_LINEAR.pdf", plot = ., width = 15, height = 25, device = "pdf")


# UNADJUSTED ####
# BOTH LINEAR AND NON-LINEAR WITH SCATTER
# this is the only plot where we exclude some outliers to make the fitted lines better visible
plot_list2 = list()
for (i in 1:(length(metrics_list))){
  temp_df = data.frame(pred = (plot_df[[i]]),age = plot_df$Age)
  #temp_df$pred = scale(temp_df$pred)
  yname=betternames_per_hemi[i]
  #names(temp_df) = c("age", "pred")
  plot_list2[[i]] =  ggplot(temp_df, aes(age,pred)) +
    geom_bin2d(bins = 150)+
    scale_fill_continuous(type = "viridis")+
    stat_smooth(method = "gam",formula = y ~ s(x, k = 4), col = "blue", se = TRUE) +
    stat_smooth(method = "lm", col = "red", se = TRUE) +
    labs(x = "Age", y = yname)+
    theme_bw()+ theme(legend.position = "none")
}
# plot the uncorrected relationships
age_curves = ggarrange(plotlist = plot_list2) %>% 
  ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_NEW_Scatter_Supplement.pdf", plot = ., width = 15, height = 25, device = "pdf") %>%
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_Scatter_Supplement.pdf", plot = ., width = 15, height = 25, device = "pdf")

#################
################# NON-LINEAR PLOTTING ####
# ADJUSTED (hemispheres stratification in one panel) ######

plot_list = list()
for (i in 1:(length(some_metrics))){
  formula1 = paste(some_metrics[i], "~ s(Age, k = 4, by = Hemisphere)+Sex+Scanner")
  yname=betternames[i]
  mod = bam(as.formula(formula1), data = both, method = "REML")
  plot_list[[i]] = plot_predictions(mod, condition = c("Age", "Hemisphere"), rug = FALSE, draw = TRUE, vcov = TRUE) + 
    theme_bw()+ labs(y = yname)
}
# plot for tsd
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_NEW_SPLINES_hemi.pdf", plot = ., width = 15, height = 25, device = "pdf")
# plot for export
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_SPLINES_hemi.pdf", plot = ., width = 15, height = 25, device = "pdf")

# ADJUSTED (sex hemisphere stratification together in one panel) ######
for (i in 1:(length(some_metrics))){
  formula1 = paste(some_metrics[i], "~ s(Age, k = 4, by = interaction(Hemisphere,Sex))+Scanner")
  yname=betternames[i]
  mod = bam(as.formula(formula1), data = both, method = "REML")
  plot_list[[i]] = plot_predictions(mod, condition = c("Age", "Hemisphere", "Sex"), rug = FALSE, draw = TRUE, vcov = TRUE) + 
    theme_bw()+ labs(y = yname)
}
# plot for tsd
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_NEW_SPLINES_sex_hemi.pdf", plot = ., width = 15, height = 25, device = "pdf")
# plot for export
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_SPLINES_sex_hemi.pdf", plot = ., width = 15, height = 25, device = "pdf")

# ADJUSTED by sex (plot each hemisphere separately) ######
# estimate gams for global scores and plot them (by sex)
plot_list = list()
for (i in 1:(length(metrics_list))){
  formula1 = paste(metrics_list[i], "~ s(Age, k = 4)+Scanner")
  yname=betternames_per_hemi[i]
  mod = bam(as.formula(formula1), data = orig_df, method = "REML")
  plot_list[[i]] = plot_predictions(mod, condition = c("Age", "Sex"), rug = FALSE, draw = TRUE, vcov = TRUE) + 
    theme_bw()+ labs(y = yname)
}

# plot for tsd
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_NEW_SPLINES_sex.pdf", plot = ., width = 15, height = 25, device = "pdf")
# plot for export
age_curves = ggarrange(plotlist = plot_list, common.legend = TRUE) %>% 
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_SPLINES_sex.pdf", plot = ., width = 15, height = 25, device = "pdf")


############### AGE-SENSITIVITY #####
# linear and non-linear age dependencies
### BOTH SEXES TOGETHER #####
# LINEAR EFFECTS ####
mod0list = list()
mod1list = list()
SS = c()
F.stat = c()
p = c()
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~",i, "+Sex+Scanner"))
  mod0list[[i]] = lm(f0,data = orig_df)
  mod1list[[i]] = lm(f1,data = orig_df)
  l = anova(mod0list[[i]], mod1list[[i]])
  SS[i] = l$`Sum of Sq`[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p = nrow(linear_hemi_effects)*p
# make df and save
linear_hemi_effects = data.frame(metrics_list, SS, F.stat, p)
linear_hemi_effects[2:4] = round(linear_hemi_effects[2:4], digits = 2)
write.csv(linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_linear_hemi_effects.csv")

# NON-LINEAR EFFECTS ####
mod3list = list()
mod4list = list()
Deviance = c()
F.stat = c()
p = c()
for (i in (metrics_list)){
  f0 = formula(paste("Age~Sex+Scanner"))
  f1 = formula(paste("Age~s(",i, ",k=4)+Sex+Scanner"))
  mod3list[[i]] = gam(f0,data = orig_df, method = "REML")
  mod4list[[i]] = gam(f1,data = orig_df, method = "REML")
  l = anova.gam(mod3list[[i]], mod4list[[i]], test = "F")
  Deviance[i] = l$Deviance[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
}
# adjust p for multiple comparison
p = nrow(non_linear_hemi_effects)*p
# make df and save
non_linear_hemi_effects = data.frame(metrics_list, Deviance, F.stat, p)
non_linear_hemi_effects[2:4] = round(non_linear_hemi_effects[2:4], digits = 2)
write.csv(non_linear_hemi_effects, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_non_linear_hemi_effects.csv")

### SEX STRATIFIED #####
males = orig_df %>% filter(Sex == "Male")
females = orig_df %>% filter(Sex == "Female")
sexes = list(males, females)
# LINEAR EFFECTS ####
mod0list = list()
mod1list = list()
SS = c()
F.stat = c()
p = c()
linear_hemi_effects_sex = list()
for (df in 1:length(sexes)){
  for (i in (metrics_list)){
  f0 = formula(paste("Age~Scanner"))
  f1 = formula(paste("Age~",i, "+Scanner"))
  mod0list[[i]] = lm(f0,data = sexes[[df]])
  mod1list[[i]] = lm(f1,data = sexes[[df]])
  l = anova(mod0list[[i]], mod1list[[i]])
  SS[i] = l$`Sum of Sq`[2]
  F.stat[i] = l$F[2]
  p[i] = l$`Pr(>F)`[2]
  }
# adjust p for multiple comparison
  p_adj = (length(p))*p
  linear_hemi_effects_sex[[df]] = data.frame(metrics_list, SS, F.stat, p, p_adj)
}
# quick check of non-age sensitive features
linear_hemi_effects_sex[[1]] %>% filter(p_adj > 0.05) #males
linear_hemi_effects_sex[[2]] %>% filter(p_adj > 0.05) #females
linear_sexes = rbindlist(linear_hemi_effects_sex)
linear_sexes$Sex = c(replicate(nrow(linear_sexes)/2, "Male"),replicate(nrow(linear_sexes)/2, "Female"))
write.csv(linear_sexes, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_linear_hemi_effects_by_sex.csv")

# NON-LINEAR EFFECTS ####
mod3list = list()
mod4list = list()
Deviance = c()
F.stat = c()
p = c()
non_linear_hemi_effects_sex = list()
for (df in 1:length(sexes)){
  for (i in (metrics_list)){
    f0 = formula(paste("Age~Scanner"))
    f1 = formula(paste("Age~s(",i, ",k=4)+Scanner"))
    mod3list[[i]] = gam(f0,data = orig_df, method = "REML")
    mod4list[[i]] = gam(f1,data = orig_df, method = "REML")
    l = anova.gam(mod3list[[i]], mod4list[[i]], test = "F")
    Deviance[i] = l$Deviance[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
  }
  # adjust p for multiple comparison
  p_adj = (length(p))*p
  non_linear_hemi_effects_sex[[df]] = data.frame(metrics_list, Deviance, F.stat, p, p_adj)
}

# quick check of non-age sensitive features
non_linear_hemi_effects_sex[[1]] %>% filter(p_adj > 0.05) #males
non_linear_hemi_effects_sex[[2]] %>% filter(p_adj > 0.05) #females
non_linear_sexes = rbindlist(non_linear_hemi_effects_sex)
non_linear_sexes$Sex = c(replicate(nrow(non_linear_sexes)/2, "Male"),replicate(nrow(non_linear_sexes)/2, "Female"))
write.csv(non_linear_sexes, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_non_linear_hemi_effects_by_sex.csv")

#########################
# HEMISPHERIC DIFFERENCES ####
# redo some of the data frames to avoid scaled frames
LEFT = raw %>% dplyr::select(contains("lh"))
RIGHT = raw %>% dplyr::select(contains("rh"))
names(LEFT)=some_metrics
names(RIGHT)=names(LEFT)

# BOTH SEXES ####
# t-tests
ttres = list()
for (i in 1:ncol(LEFT)){
  x = LEFT[i]
  y = RIGHT[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  ttres[[i]] = t.test(df$x, df$y, paired = T)
}
ttres_p = c()
name = c()
for (i in 1:length(ttres)){
  ttres_p[i] = (ttres[[i]]$p.value)
}
tdat = data.frame(some_metrics, ttres_p)
# absolute value
sum(ifelse(ttres_p > 0.05/length(ttres), 1, 0))
# relative value
sum(ifelse(ttres_p > 0.05/length(ttres), 1, 0))/length(ttres)
# show which values are non-sig
tdat %>% filter(ttres_p > 0.05/length(ttres))

# SEX-STRATIFIED ####
LEFT_m = raw %>% dplyr::select(contains("lh"), Sex) %>% filter(Sex == "Male") %>% dplyr::select(!Sex)
RIGHT_m = raw %>% dplyr::select(contains("rh"), Sex) %>% filter(Sex == "Male") %>% dplyr::select(!Sex)
names(LEFT_m)=some_metrics
names(RIGHT_m)=some_metrics
LEFT_f = raw %>% dplyr::select(contains("lh"), Sex) %>% filter(Sex == "Female") %>% dplyr::select(!Sex)
RIGHT_f = raw %>% dplyr::select(contains("rh"), Sex) %>% filter(Sex == "Female") %>% dplyr::select(!Sex)
names(LEFT_f)=some_metrics
names(RIGHT_f)=some_metrics

ttres_m = list()
ttres_f = list()
for (i in 1:ncol(LEFT_m)){
  x = LEFT_m[i]
  y = RIGHT_m[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  ttres_m[[i]] = t.test(df$x, df$y, paired = T)
  x = LEFT_f[i]
  y = RIGHT_f[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  ttres_f[[i]] = t.test(df$x, df$y, paired = T)
}
ttres_p_m = c()
ttres_p_f = c()

name = c()
for (i in 1:length(ttres_m)){
  ttres_p_m[i] = (ttres_m[[i]]$p.value)
  ttres_p_f[i] = (ttres_f[[i]]$p.value)
}
tdat_m = data.frame(some_metrics, ttres_p_m)
tdat_f = data.frame(some_metrics, ttres_p_f)
#
#
# SHOW RESULTS
#
# for males
#
# absolute value
sum(ifelse(ttres_p_m > 0.05/length(ttres_m), 1, 0))
# relative value
sum(ifelse(ttres_p_m > 0.05/length(ttres_m), 1, 0))/length(ttres)
# show which values are non-sig
tdat_m %>% filter(ttres_p_m > 0.05/length(ttres_m))
#
# for females
#
# absolute value
sum(ifelse(ttres_p_f > 0.05/length(ttres_f), 1, 0))
# relative value
sum(ifelse(ttres_p_f > 0.05/length(ttres_f), 1, 0))/length(ttres)
# show which values are non-sig
tdat_f %>% filter(ttres_p_f > 0.05/length(ttres_f))


