# This script estimates AGE SENSITIVITY of ABSOLUTE brain asymmetries BY SEX
# It also compares linear and non-linear models using AIC and BIC
# Date: 02.11.2023 
# Author: Max Korbmacher (max.korbmacher@gmail.com)
# 
print("This script estimates AGE SENSITIVITY of ABSOLUTE regional brain asymmetries BY SEX.")
print("It also compares linear and non-linear models using AIC and BIC")
print("Date: 02.11.2023 ")
print("Author: Max Korbmacher (max.korbmacher@gmail.com)")


# PREP ####
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, marginaleffects, ggpubr, ggplot2, Rmpfr, mgcv, lsr)
# load data
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/T1.csv")
dMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/dMRI.csv")
multimodal = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/combi.csv")
demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")

# split left from right
LT1 = T1 %>% dplyr::select(eid, starts_with(c("lh", "Left")))
RT1 = T1 %>% dplyr::select(eid, starts_with(c("rh","Right")))
LdMRI = dMRI %>% dplyr::select(eid, ends_with(c("L")))
RdMRI = dMRI %>% dplyr::select(eid, ends_with(c("R")))
L = merge(LT1, LdMRI, by = "eid")
R = merge(RT1,RdMRI, by = "eid")

# add demo to L for later regression
L = merge(L,demo,by="eid")
LdMRI = merge(LdMRI, demo, by = "eid")
LT1 = merge(LT1,demo,by="eid")

## for multimodal data
L = L %>% dplyr::select(-eid)
R = R %>% dplyr::select(-eid)

## dMRI only
LdMRI = LdMRI %>% dplyr::select(-eid)
RdMRI = RdMRI %>% dplyr::select(-eid)
## T1 only
LT1 = LT1 %>% dplyr::select(-eid)
RT1 = RT1 %>% dplyr::select(-eid)

# define laterality index/asymmetry function (ABSOLUTE VALUES)
LI = function(L,R){
  abs((L-R)/(L+R))
}
# now, produce mean-centered (to zero), scaled laterality indexed regional metrics
LIdat = data.frame(matrix(ncol = ncol(R),nrow = nrow(R)))
LIdMRI = data.frame(matrix(ncol = ncol(RdMRI),nrow = nrow(RdMRI)))
LIT1 = data.frame(matrix(ncol = ncol(RT1),nrow = nrow(RT1)))
for (i in 1:ncol(R)){
  LIdat[i] = as.numeric(scale(LI(L[[i]], R[[i]])))
}
for (i in 1:ncol(RdMRI)){
  LIdMRI[i] = as.numeric(scale(LI(LdMRI[i], RdMRI[i])))
}
for (i in 1:ncol(RT1)){
  LIT1[i] = as.numeric(scale(LI(LT1[i], RT1[i])))
}
names(LIdat) = names(R)
names(LIdMRI) = names(RdMRI)
names(LIT1) = names(RT1)

# remove mean average remainders
remainders = function(data){
  data %>% dplyr::select(!contains("mean"), !contains("total"))
}
LIdat = remainders(LIdat)
LIdMRI = remainders(LIdMRI)
LIT1 = remainders(LIT1)
#
#
# add sex, age, site
LIdat$sex = L$Sex
LIdat$age = L$Age
LIdat$site = L$Scanner
LIdMRI$sex = LdMRI$Sex
LIdMRI$age = LdMRI$Age
LIdMRI$site = LdMRI$Scanner
LIT1$sex = LT1$Sex
LIT1$age = LT1$Age
LIT1$site = LT1$Scanner
#
# scale age for easier plotting of coefficients
#LIdat$age_scaled = scale(LIdat$age)
LIdMRI$age_scaled = scale(LIdMRI$age)
LIT1$age_scaled = scale(LIT1$age)
# split data frames by sex
LIdMRI_m = LIdMRI %>% filter(sex == "Male")
LIdMRI_f = LIdMRI %>% filter(sex == "Female")
LIT1_m = LIT1 %>% filter(sex == "Male")
LIT1_f = LIT1 %>% filter(sex == "Female")
#
# make metrics list containing only MRI variables
a = LIdMRI_m %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
b = LIdMRI_f %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
c = LIT1_m %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
d = LIT1_f %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
metrics_list = list(a,b,c,d)
#
# make a list of the relevant data frames
LRTlist = list(LIdMRI_m, LIdMRI_f, LIT1_m, LIT1_f)
#
# TEST AGE SENSITIVITY BOTH SEXES #####
#
### linear modeling #### 
linear_hemi_effects = list()
print("Starting linear modelling.")
for (df in 1:length(LRTlist)){
  dat = LRTlist[[df]]
  dat = na.omit(dat)
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  AICf = c()
  BICf = c()
  for (i in (metrics_list[[df]])){
    f0 = formula(paste("age_scaled~site"))
    f1 = formula(paste("age_scaled~",i, "+site"))
    mod0list = lm(f0,data = dat)
    mod1list = lm(f1,data = dat)
    l = anova(mod0list, mod1list)
    SS[i] = l$`Sum of Sq`[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    AICf[i] = AIC(mod1list)
    BICf[i] = BIC(mod1list)
    print(paste("LRT assessing the LINEAR age-sensitivity of LI of ", i, "completed."))
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], SS, F.stat, p, p_adj, AICf, BICf)
  print(paste("Data set ", df, " processed."))
}
print("Writing files")
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_dMRI_MALES.csv")
print("...")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_dMRI_FEMALES.csv")
print("...")
write.csv(linear_hemi_effects[[3]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_T1_MALES.csv")
print("...")
write.csv(linear_hemi_effects[[4]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_T1_FEMALES.csv")
print("Writing files done")
#
print("Starting non-linear modelling.")
### non-linear modeling ####
non_linear_hemi_effects = list()
for (df in 1:length(LRTlist)){
  dat = LRTlist[[df]]
  dat = na.omit(dat)
  SS = c()
  F.stat = c()
  p = c()
  Deviance = c()
  p_adj = c()
  AICf = c()
  BICf = c()
  for (i in (metrics_list[[df]])){
    f0 = formula(paste("age_scaled~site"))
    f1 = formula(paste("age_scaled~s(",i, ",k=4)+site"))
    mod0list = gam(f0,data = dat, method = "REML")
    mod1list = gam(f1,data = dat, method = "REML")
    l = anova.gam(mod0list, mod1list, test = "F")
    Deviance[i] = l$Deviance[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    AICf[i] = AIC(mod1list)
    BICf[i] = BIC(mod1list)
    print(paste("LRT assessing the age-sensitivity of NON-LINEAR LI in", i, "completed."))
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  non_linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], Deviance, F.stat, p, p_adj, AICf, BICf)
  paste("Data set ", df, " processed.")
}
print("Writing files")
write.csv(non_linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_dMRI_MALES.csv")
print("...")
write.csv(non_linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_dMRI_FEMALES.csv")
print("...")
write.csv(non_linear_hemi_effects[[3]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_T1_MALES.csv")
print("...")
write.csv(non_linear_hemi_effects[[4]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_T1_FEMALES.csv")
print("Writing files done.")
#
#
#
#
#
print("####")
print("####")
print("###################################")
print("###################################")
print("Analyses completed for each males and females.")
print("###################################")
print("###################################")
print("The following assesses the distribution of AIC and BIC of the full linear vs non-linear models of age ~ feature + sex + scanner")
print("NOTE: p-values need to be adjusted as multiple tests are executed on the same hypothesis!")
print("#")
print("Differences in model fit for diffusion-weighted features IN MALES:")
print("###################################")
print("Start with AIC paired t-test:")
t = t.test(linear_hemi_effects[[1]]$AICf, non_linear_hemi_effects[[1]]$AICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[1]]$AICf, non_linear_hemi_effects[[1]]$AICf, method = "paired")
print(d)
print("###################################")
print("BIC paired t-test:")
t = t.test(linear_hemi_effects[[1]]$BICf, non_linear_hemi_effects[[1]]$BICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[1]]$BICf, non_linear_hemi_effects[[1]]$BICf, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("####")
print("Differences in model fit for diffusion-weighted features IN FEMALE:")
print("###################################")
print("Start with AIC paired t-test:")
t = t.test(linear_hemi_effects[[2]]$AICf, non_linear_hemi_effects[[2]]$AICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[2]]$AICf, non_linear_hemi_effects[[2]]$AICf, method = "paired")
print(d)
print("###################################")
print("BIC paired t-test:")
t = t.test(linear_hemi_effects[[2]]$BICf, non_linear_hemi_effects[[2]]$BICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[2]]$BICf, non_linear_hemi_effects[[2]]$BICf, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("####")
print("###################################")
print("Differences in model fit for T1-weighted features IN MALES:")
print("###################################")
print("Start with AIC paired t-test:")
t = t.test(linear_hemi_effects[[3]]$AICf, non_linear_hemi_effects[[3]]$AICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[3]]$AICf, non_linear_hemi_effects[[3]]$AICf, method = "paired")
print(d)
print("###################################")
print("BIC paired t-test:")
t = t.test(linear_hemi_effects[[3]]$BICf, non_linear_hemi_effects[[3]]$BICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[3]]$BICf, non_linear_hemi_effects[[3]]$BICf, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("####")
print("####")
print("####")
print("###################################")
print("Differences in model fit for T1-weighted features IN FEMALES:")
print("###################################")
print("Start with AIC paired t-test:")
t = t.test(linear_hemi_effects[[4]]$AICf, non_linear_hemi_effects[[4]]$AICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[4]]$AICf, non_linear_hemi_effects[[4]]$AICf, method = "paired")
print(d)
print("###################################")
print("BIC paired t-test:")
t = t.test(linear_hemi_effects[[4]]$BICf, non_linear_hemi_effects[[4]]$BICf, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[4]]$BICf, non_linear_hemi_effects[[4]]$BICf, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("Now, we present percentages of significantly age-associated LI features")
print("####")
ress = c("dMRI lin males", "dMRI lin females", "dMRI NON-lin males", "dMRI NON-lin females", "t1 lin males", "t1 lin females", "t1 NON-lin males", "t1 NON-lin females")
print("Starting with NON-LINEAR EFFECTS")
for (i in 1:length(non_linear_hemi_effects)){
  print(ress[i])
  p_adj = non_linear_hemi_effects[[i]]$p_adj
  print("Total number of non-sig/non-age sens features:")
  print(paste(nrow(non_linear_hemi_effects[[i]])-sum(ifelse(p_adj > 0.05/nrow(non_linear_hemi_effects[[i]]), 1, 0)),"of", nrow(non_linear_hemi_effects[[i]]), "were significantly age-related using non-linear associations"))
  #print("Percentages:")
  print(paste("Proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/nrow(non_linear_hemi_effects[[i]]), 1, 0))/nrow(non_linear_hemi_effects[[i]]), "percent."))
  print("mean absolute F-statistic (SD)")
  print(non_linear_hemi_effects[[i]] %>% filter(p_adj < 0.05/nrow(non_linear_hemi_effects[[i]])) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat))))
  print("#########################")
}
print("####")
print("NOW, LINEAR EFFECTS")
for (i in 1:length(linear_hemi_effects)){
  print(ress[i])
  p_adj = linear_hemi_effects[[i]]$p_adj
  print("Total number of non-sig/non-age sens features:")
  print(paste(nrow(linear_hemi_effects[[i]])-sum(ifelse(p_adj > 0.05/nrow(linear_hemi_effects[[i]]), 1, 0)),"of", nrow(linear_hemi_effects[[i]]), "were significantly age-related using non-linear associations"))
  #print("Percentages:")
  print(paste("Proportion of significantly age sensitive metrics = ", 1-sum(ifelse(p_adj > 0.05/nrow(linear_hemi_effects[[i]]), 1, 0))/nrow(linear_hemi_effects[[i]]), "percent."))
  print("mean absolute F-statistic (SD)")
  print(linear_hemi_effects[[i]] %>% filter(p_adj < 0.05/nrow(linear_hemi_effects[[i]])) %>% summarize(Fmean = mean(abs(F.stat)), Fsd = sd(abs(F.stat))))
  print("#########################")
}
print("End of script. Job done.")