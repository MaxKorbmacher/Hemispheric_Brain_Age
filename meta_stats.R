# meta-stats for 
## a) regional L-R hemisphere differences in T1w and diffusion metrics
## b) age sensitivity

# load packages
library(dplyr)
library(ggplot2)
library(ggpubr)

# load data
#t1_diff = read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_T1w_features_diff.csv")
#dmri_diff = read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_dMRI_features_diff.csv")
#t1_age = read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_T1w_features_age_sens.csv")
#dmri_age = read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_dMRI_features_agesens.csv")
setwd("/Users/max/Documents/Projects/Hemi_brain_age/")
t1_LI = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_both_sexes.csv")
dmri_LI =  read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_both_sexes.csv")
t1_LI_m = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_MALES.csv")
t1_LI_f= read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_FEMALES.csv")
dmri_LI_m = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_MALES.csv")
dmri_LI_f = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_FEMALES.csv")
# # show distributions of effects
# ############# 1) for feature asymmetries
# ## T1w MRI
# t1_diff %>% filter(T1_p_adj < .05) %>% summarise(Mean = mean(T1_CohensD), SD = sd(T1_CohensD))
# ## diffusion MRI
# dmri_diff %>% filter(dMRI_p_adj < .05) %>% summarise(Mean = mean(dMRI_CohensD), SD = sd(dMRI_CohensD))
# ########### 2) for features' age-dependencies
# ## T1w MRI
# t1_age %>% filter(ps < .05) %>% summarise(Mean = mean(F_), SD = sd(F_))
# ## diffusion MRI
# dmri_age %>% filter(ps < .05) %>% summarise(Mean = mean(F_), SD = sd(F_))
########### 3) for LI age-dependencies
# estimate adjusted p-vals
t1_LI$p_adj = (2 * pnorm(-abs(t1_LI$t)))*nrow(t1_LI)
t1_LI_m$p_adj = (2 * pnorm(-abs(t1_LI_m$t)))*nrow(t1_LI_m)
t1_LI_f$p_adj = (2 * pnorm(-abs(t1_LI_f$t)))*nrow(t1_LI_f)
dmri_LI$p_adj = (2 * pnorm(-abs(dmri_LI$t)))*nrow(dmri_LI)
dmri_LI_m$p_adj = (2 * pnorm(-abs(dmri_LI_m$t)))*nrow(dmri_LI_m)
dmri_LI_f$p_adj = (2 * pnorm(-abs(dmri_LI_f$t)))*nrow(dmri_LI_f)

# show mean and SD of the effect sizes
t1_LI %>% filter(t1_LI$p_adj < .05) %>% summarise(Mean = mean(abs(beta)), SD = sd(abs(beta)))
dmri_LI %>% filter(dmri_LI$p_adj < .05) %>% summarise(Mean = mean(abs(beta)), SD = sd(abs(beta)))
#
# # focus on mean values only
# betas = c(0.026,.01,.068,0.018,.0250,.086,.102,.054,.022,.037,.022,.026,.027,.03,.052,.099,.083,.072,.04,.097,.056)
# mean(betas)
# sd(betas)
#
# Density plots for all slopes (incl non sig.) ####
# dMRI
dmri_both = ggplot(dmri_LI, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Both sexes' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
dmri_m = ggplot(dmri_LI_m, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Males' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
dmri_f = ggplot(dmri_LI_f, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Females' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
# T1w
t1_both = ggplot(t1_LI, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Both sexes' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
t1_m = ggplot(t1_LI_m, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Males' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
t1_f = ggplot(t1_LI_f, aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Females' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
uncorrected = ggarrange(dmri_both, dmri_m, dmri_f, t1_both, t1_m, t1_f)
ggsave('all_slopes.pdf', uncorrected, width = 17, height = 9 ) 

# Density plots for all slopes (incl ONLY sig.) ####
# dMRI
dmri_both = dmri_LI %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Both sexes' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()

dmri_m = dmri_LI_m %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Males' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
dmri_f = dmri_LI_f %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Females' dMRI features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
# T1w
t1_both = t1_LI %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) + 
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Both sexes' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
t1_m = t1_LI_m %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) + 
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Males' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
t1_f = t1_LI_f %>% filter(p_adj < .05) %>% ggplot(aes(x=beta)) + 
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",
             linetype="dashed")+
  labs(title="Females' T1 features' age associations",x="Slopes", y = "Density")+
  geom_boxplot() +   theme_classic()
corrected = ggarrange(dmri_both, dmri_m, dmri_f, t1_both, t1_m, t1_f)
ggsave('sig_slopes.pdf', corrected, width = 17, height = 9 ) 
