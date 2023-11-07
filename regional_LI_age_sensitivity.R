# This script estimates AGE SENSITIVITY of ABSOLUTE brain asymmetries
# It also compares linear and non-linear models using AIC and BIC
# Date: 30.10.2023 
# Author: Max Korbmacher (max.korbmacher@gmail.com)
# 
print("This script estimates AGE SENSITIVITY of ABSOLUTE regional brain asymmetries.")
print("It also compares linear and non-linear models using AIC and BIC")
print("Date: 30.10.2023 ")
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
# make metrics list containing only MRI variables
metrics_list = list(names(LIdMRI), names(LIT1))
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
LIdat$age_scaled = scale(LIdat$age)
LIdMRI$age_scaled = scale(LIdMRI$age)
LIT1$age_scaled = scale(LIT1$age)
# make list of LI data frames
LIdfs = list(LIdat, LIdMRI, LIT1)
#
# TEST AGE SENSITIVITY BOTH SEXES #####
# make a list of the relevant data frames
LRTlist = list(LIdMRI, LIT1)
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
  aic = c()
  bic = c()
  for (i in (metrics_list[[df]])){
    f0 = formula(paste("age_scaled~site"))
    f1 = formula(paste("age_scaled~",i, "+site"))
    mod0list = lm(f0,data = dat)
    mod1list = lm(f1,data = dat)
    l = anova(mod0list, mod1list)
    SS[i] = l$`Sum of Sq`[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    aic[i] = AIC(mod1list)
    bic[i] = BIC(mod1list)
    print(paste("LRT assessing the LINEAR age-sensitivity of", i, "completed."))
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], SS, F.stat, p, p_adj, aic, bic)
  print(paste("Data set ", df, " processed."))
}
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_dMRI.csv")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_sens_T1.csv")
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
  aic = c()
  bic = c()
  for (i in (metrics_list[[df]])){
    f0 = formula(paste("age_scaled~site"))
    f1 = formula(paste("age_scaled~s(",i, ",k=4)+site"))
    mod0list = gam(f0,data = dat, method = "REML")
    mod1list = gam(f1,data = dat, method = "REML")
    l = anova.gam(mod0list, mod1list, test = "F")
    Deviance[i] = l$Deviance[2]
    F.stat[i] = l$F[2]
    p[i] = l$`Pr(>F)`[2]
    aic[i] = AIC(mod1list)
    bic[i] = BIC(mod1list)
    paste("LRT assessing the age-sensitivity of NON-LINEAR", i, "completed.")
  }
  # adjust p for multiple comparison
  p_adj = length(p)*p
  # make df and save
  non_linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], Deviance, F.stat, p, p_adj, aic, bic)
  paste("Data set ", df, " processed.")
}
write.csv(non_linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_dMRI.csv")
write.csv(non_linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_NON_LINEAR_LI_age_sens_T1.csv")
#
print("####")
print("####")
print("###################################")
print("###################################")
print("Analyses completed for both sexes.")
print("###################################")
print("###################################")
print("The following assesses the distribution of AIC and bic of the full linear vs non-linear models of age ~ feature + sex + scanner")
print("#")
print("Differences in model fit for diffusion-weighted features")
print("###################################")
print("Start with AIC paired t-test:")
t = t.test(linear_hemi_effects[[1]]$aic, non_linear_hemi_effects[[1]]$aic, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[1]]$aic, non_linear_hemi_effects[[1]]$aic, method = "paired")
print(d)
print("###################################")
print("bic paired t-test:")
t = t.test(linear_hemi_effects[[1]]$bic, non_linear_hemi_effects[[1]]$bic, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[1]]$bic, non_linear_hemi_effects[[1]]$bic, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("####")
print("###################################")
print("Differences in model fit for T1-weighted features:")
print("###################################")
print("Start with aic paired t-test:")
t = t.test(linear_hemi_effects[[2]]$aic, non_linear_hemi_effects[[2]]$aic, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[2]]$aic, non_linear_hemi_effects[[2]]$aic, method = "paired")
print(d)
print("###################################")
print("bic paired t-test:")
t = t.test(linear_hemi_effects[[2]]$bic, non_linear_hemi_effects[[2]]$bic, method = "paired")
print(t)
print("####")
print("Effect size, assuming paired/dependent samples:")
d = cohensD(linear_hemi_effects[[2]]$bic, non_linear_hemi_effects[[2]]$bic, method = "paired")
print(d)
print("###################################")
print("####")
print("####")
print("We can also calucalate how many features were age-sensitive.")
print("###################################")
print("Age sensitivity of diffusion-weighted features:")
print("###################################")
print(paste("first non-linear then linear model's estimation of age-sensitivity"))
print(non_linear_hemi_effects[[1]] %>% filter(p_adj < .05) %>% summarize(N_sig = length(p_adj), N_rel = length(p_adj)/nrow(non_linear_hemi_effects[[1]])))
print(linear_hemi_effects[[1]] %>% filter(p_adj < .05) %>% summarize(N_sig = length(p_adj), N_rel = length(p_adj)/nrow(non_linear_hemi_effects[[1]])))
print("###################################")
print("Age sensitivity of T1-weighted features:")
print("###################################")
print(non_linear_hemi_effects[[2]] %>% filter(p_adj < .05) %>% summarize(N_sig = length(p_adj), N_rel = length(p_adj)/nrow(non_linear_hemi_effects[[1]])))
print(linear_hemi_effects[[2]] %>% filter(p_adj < .05) %>% summarize(N_sig = length(p_adj), N_rel = length(p_adj)/nrow(non_linear_hemi_effects[[1]])))

print("End of script. Job done.")