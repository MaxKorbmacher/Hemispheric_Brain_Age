# This script estimates slopes for the (linear and non-linear) relationship between brain asymmetries and age
# created 05.09.2023, updated 1.11.2023, Max Korbmacher (max.korbmacher@gmail.com)

# load packages
library(dplyr)
library(ggplot2)
library(ggpubr)

# load data
L = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Lcombi.csv")
R = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Rcombi.csv")
LT1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/LT1.csv")
RT1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/RT1.csv")
LdMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/LdMRI.csv")
RdMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/RdMRI.csv")
demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")

# add demo to L for later regression
L = merge(L,demo,by="eid")
LdMRI = merge(LdMRI, demo, by = "eid")
LT1 = merge(LT1,demo,by="eid")

# remove counter and eid
## for multimodal data
L = L %>% select(-X,-eid)
R = R %>% select(-X,-eid)
## dMRI only
LdMRI = LdMRI %>% select(-X,-eid, age)
RdMRI = RdMRI %>% select(-X,-eid, age)
## T1 only
LT1 = LT1 %>% select(-X,-eid, age)
RT1 = RT1 %>% select(-X,-eid, age)

# define laterality index/asymmetry function (ABSOLUTE VALUES)
LI = function(L,R){
  abs((L-R)/(L+R))
}

# now, produce mean-centered (to zero), scaled laterality indexed regional metrics
LIdat = data.frame(matrix(ncol = ncol(R),nrow = nrow(R)))
LIdMRI = data.frame(matrix(ncol = ncol(RdMRI),nrow = nrow(RdMRI)))
LIT1 = data.frame(matrix(ncol = ncol(RT1),nrow = nrow(RT1)))
for (i in 1:ncol(R)){
  LIdat[i] = scale(LI(L[i], R[i]))
}
for (i in 1:ncol(RdMRI)){
  LIdMRI[i] = scale(LI(LdMRI[i], RdMRI[i]))
}
for (i in 1:ncol(RT1)){
  LIT1[i] = scale(LI(LT1[i], RT1[i]))
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

# scale age for easier plotting of coefficients
LIdat$age_scaled = scale(LIdat$age)
LIdMRI$age_scaled = scale(LIdMRI$age)
LIT1$age_scaled = scale(LIT1$age)

# test age associations
# (due to convergence problems, we use simple linear and generalized additive models, instead of site as random effect)
## multimodal
coeff = c()
pvals = c()
for (i in 1:ncol(R)){
  f = formula(paste(colnames(LIdat)[i],"~age_scaled+sex+site"))
  tmpmodel = lm(f, data = LIdat)
  coeff[i] = summary(tmpmodel)$coefficients[2]
  pvals[i] = summary(tmpmodel)$coefficients[2,4]
}
regional = data.frame(names(R),coeff,pvals)
## diffusion
coeff = c()
pvals = c()
for (i in 1:ncol(RdMRI)){
  f = formula(paste(colnames(LIdMRI)[i],"~age_scaled+sex+site"))
  tmpmodel = lm(f, data = LIdMRI)
  coeff[i] = summary(tmpmodel)$coefficients[2]
  pvals[i] = summary(tmpmodel)$coefficients[2,4]
}
regional_dMRI = data.frame(names(RdMRI),coeff,pvals)
## T1
coeff = c()
pvals = c()
for (i in 1:ncol(RT1)){
  f = formula(paste(colnames(LIT1)[i],"~age_scaled+sex+site"))
  tmpmodel = lm(f, data = LIT1)
  coeff[i] = summary(tmpmodel)$coefficients[2]
  pvals[i] = summary(tmpmodel)$coefficients[2,4]
}
regional_T1 = data.frame(names(RT1),coeff,pvals)


## assess the results ####

####### OVERALL / MULTIMODAL

# of 958 asymmetries
nrow(regional)
# 582 are significantly related with age after Bonferroni correction
plot_reg = regional %>% dplyr::filter(pvals < 0.05/nrow(regional))
nrow(plot_reg)
#
#
# !!! WATCH OUT !!!
# for plotting, I removed one extreme finding in the accumbens (strong asymmetry increase)
plot_reg = plot_reg %>% dplyr::filter(!names.R. == "Right.Accumbens.area")

# density plot of the effects 
plot1 = ggplot(plot_reg, aes(x=coeff)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="Multimodal metrics",x="Asymmetry-age relationship", y = "Density")+
  theme_classic()
# check meta stats
regional %>% dplyr::filter(pvals < 0.05/nrow(regional)) %>% summarise(M = mean(coeff), SD = sd(coeff), MD = median(coeff), MAD = mad(coeff))

####### DIFFUSION
regional_dMRI = regional_dMRI %>% dplyr::filter(!names.RdMRI. == "age")
# of 840 asymmetries
nrow(regional_dMRI)
# 436 are significantly related with age after Bonferroni correction
plot_dMRI = regional_dMRI %>% dplyr::filter(pvals < 0.05/nrow(regional_dMRI))
nrow(plot_dMRI)
# density plot of the effects 
plot2 = ggplot(plot_dMRI, aes(x=coeff)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="Diffusion Metrics",x="Asymmetry-age relationship", y = "Density")+
  theme_classic()
# check meta stats
regional_dMRI %>% dplyr::filter(pvals < 0.05/nrow(regional_dMRI)) %>% summarise(M = mean(coeff), SD = sd(coeff), MD = median(coeff), MAD = mad(coeff))

### T1-WEIGHTED
regional_T1 = regional_T1 %>% dplyr::filter(!names.RT1. == "age")
# of 117 asymmetries
nrow(regional_T1)
# 53 are significantly related with age after Bonferroni correction
plot_T1 = regional_T1 %>% dplyr::filter(pvals < 0.05/nrow(regional_T1))
nrow(plot_T1)
# density plot of the effects 
plot3 = ggplot(plot_T1, aes(x=coeff)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="T1 Metrics",x="Asymmetry-age relationship", y = "Density")+
  theme_classic()
# check meta stats
regional_T1 %>% dplyr::filter(pvals < 0.05/nrow(regional_dMRI)) %>% summarise(M = mean(coeff), SD = sd(coeff), MD = median(coeff), MAD = mad(coeff))

# for replication purposes:
regional_T1 %>% filter(grepl("thickness",names.RT1.))
# previous study shows decreasing asymmetry in cortical thickness, we cannot confirm this in cross-sectional data
#
#
#
# write the three data frames indicating absolute LI associations with age
write.csv(regional,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_absolute_LI_multimodal.csv")
write.csv(regional_dMRI,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_absolute_LI_dMRI.csv")
write.csv(regional_T1,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_absolute_LI_T1.csv")

#### COMBINE PLOTS
dens_plot_lin = ggarrange(plot1,0, plot2, plot3)
dens_plot_non_lin = ggarrange(plot1,0, plot2, plot3)
dens_plot = ggarrange(dens_plot_lin, dens_plot_non_lin, labels = c("a", "b"))


#ggsave("/tsd/p33/home/p33-maxk/export/Hemi_dens.pdf",dens_plot, height =8, width = 8)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_Supplement10.pdf",dens_plot, height =8, width = 8)
