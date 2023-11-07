print("This script estimates AGE ASSOCIATIONS of ABSOLUTE regional brain asymmetries BY SEX.")
print("Date: 06.11.2023 ")
print("Author: Max Korbmacher (max.korbmacher@gmail.com)")


# PREP ####
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, marginaleffects, ggpubr, ggplot2, Rmpfr, mgcv, lsr, lm.beta)
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
# make also another list for dfs where both data frames are included 
LRTlist2 = list(LIdMRI, LIT1)
#
# TEST AGE ASSOCIATIONS BOTH SEXES #####
#
### linear modeling #### 
linear_hemi_effects = list()
print("Starting linear modelling.")
for (df in 1:length(LRTlist)){
  dat = LRTlist[[df]]
  dat = na.omit(dat)
  t = c()
  beta = c()
  for (i in (metrics_list[[df]])){
    f1 = formula(paste("age_scaled~",i, "+site"))
    mod1list = lm(f1,data = dat)
    beta[i] = lm.beta(mod1list)$standardized.coefficients[2]
    t[i] = summary(mod1list)$coefficients[2,3]
    print(paste("LRT assessing the LINEAR age-association of LI of ", i, "completed."))
  }
  # adjust p for multiple comparison
  #p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], beta, t)
  print(paste("Data set ", df, " processed."))
}
print("Writing files")
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_MALES.csv")
print("file 1 done ...")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_FEMALES.csv")
print("file 2 done ...")
write.csv(linear_hemi_effects[[3]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_MALES.csv")
print("file 3 done ...")
write.csv(linear_hemi_effects[[4]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_FEMALES.csv")
print("Writing files done")
print("###")
print("Print some meta stats.")
# 
# # check meta stats
# malesDMRI = read.csv( "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_MALES.csv")
# femalesDMRI = read.csv( "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_FEMALES.csv")
# malesT1 = read.csv( "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_MALES.csv")
# femalesT1 = read.csv( "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_FEMALES.csv")
# linear_hemi_effects = list(malesDMRI, femalesDMRI, malesT1, femalesT1)
#.N <- function(.) mpfr(., precBits = 200)
# # 2 * pnorm(-abs(.N(betalist[[i]]$t)))
print("now print unadjusted standardized beta values (males dMRI, females dMRI, males T1, females T1)")
for (i in 1:4){
  print(linear_hemi_effects[[i]] %>% summarize(mean(beta), sd(beta)))
  print("######################")
}
print("now print adjusted standardized beta values (males dMRI, females dMRI, males T1, females T1)")
for (i in 1:4){
  linear_hemi_effects[[i]]$p_val = 2 * pnorm(-abs(linear_hemi_effects[[i]]$t))*nrow(linear_hemi_effects[[i]])
  print(linear_hemi_effects[[i]] %>%  filter(p_val < 0.05) %>% summarize(mean(beta), sd(beta)))
  print("######################")
} 
#
print("######################")
print("######################")
print("we also estimate the same for data of both sexes")
linear_hemi_effects = list()
print("Starting linear modelling.")
a = LIdMRI %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
b = LIT1 %>% dplyr::select(-sex, -age, -site, -age_scaled) %>% names()
metrics_list = list(a, b)
for (df in 1:length(LRTlist2)){
  dat = LRTlist2[[df]]
  dat = na.omit(dat)
  t = c()
  beta = c()
  for (i in (metrics_list[[df]])){
    f1 = formula(paste("age_scaled~",i, "+site"))
    mod1list = lm(f1,data = dat)
    beta[i] = lm.beta(mod1list)$standardized.coefficients[2]
    t[i] = summary(mod1list)$coefficients[2,3]
    print(paste("LRT assessing the LINEAR age-association of LI of ", i, "completed."))
  }
  # adjust p for multiple comparison
  #p_adj = length(p)*p
  # make df and save
  linear_hemi_effects[[df]] = data.frame(metrics_list[[df]], beta, t)
  print(paste("Data set ", df, " processed."))
}
print("Writing files")
write.csv(linear_hemi_effects[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_both_sexes.csv")
print("file 1 done ...")
write.csv(linear_hemi_effects[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_both_sexes.csv")
print("file 2 done ...")
print("now print unadjusted standardized beta values (dMRI, then T1 for both sexes)")
for (i in 1:2){
  print(linear_hemi_effects[[i]] %>% summarize(mean(beta), sd(beta)))
  print("######################")
}
print("now print adjusted standardized beta values (dMRI, then T1 for both sexes)")
for (i in 1:2){
  linear_hemi_effects[[i]]$p_val = 2 * pnorm(-abs(linear_hemi_effects[[i]]$t))*nrow(linear_hemi_effects[[i]])
  print(linear_hemi_effects[[i]] %>%  filter(p_val < 0.05) %>% summarize(mean(beta), sd(beta)))
  print("######################")
} 
print("This concludes the job.")

# check also grant mean
a = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_both_sexes.csv")
b = read.csv("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_both_sexes.csv")          
c = rbind(a,b)
c$p_val = 2 * pnorm(-abs(c$t))*nrow(c)
print(c %>%  filter(p_val < 0.05) %>% summarize(mean(beta), sd(beta)))