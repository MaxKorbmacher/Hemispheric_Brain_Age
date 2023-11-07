# exploring hemispheric age predictions
# updated: 02.11.2023 adding sex stratification

## PREP ####

# packages
#install.packages("afex")
#install.packages("tidyquant")
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
#install.packages("ppcor")
library(ppcor)
library(MuMIn)
#install.packages("ggdist")
library(ggdist)
library(emmeans)

# read in data
LT1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LT1.csv")
RT1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RT1.csv")
LdMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LdMRI.csv")
RdMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RdMRI.csv")
Lcombi = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_Lcombi.csv")
Rcombi = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_Rcombi.csv")
demo = read.csv("/cluster/projects/p33/users/maxk/UKB/brainage_MK/BAG_associations/final_dwMRI_data.csv")

bothT1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_T1.csv")
bothdMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_dwMRI.csv")
bothcombi = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_combi.csv")
# variables of interest
handedness = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/analyses/handedness.csv")
vois = read.csv("/cluster/projects/p33/users/maxk/UKB/BP_BAG/data/demo_clean.csv")


# put predictions from different hemispheres into one data set
T1_pred = data.frame(LT1$eid,LT1$age, LT1$pred_age_LT1, RT1$pred_age_RT1)
T1all = data.frame(bothT1$eid,bothT1$pred_age_T1)
dMRI_pred = data.frame(LdMRI$eid, LdMRI$pred_age_LdMRI,RdMRI$pred_age_RdMRI)
dMRIall = data.frame(bothdMRI$eid,bothdMRI$pred_age_dwMRI)
combi_pred = data.frame(Lcombi$eid, Lcombi$pred_age_Lcombi,Rcombi$pred_age_Rcombi)
combiall = data.frame(bothcombi$eid,bothcombi$pred_age_combi)


# change column names
colnames(T1_pred) = c("eid","age", "LT1", "RT1")
colnames(T1all) = c("eid", "bothT1")
colnames(dMRI_pred) = c("eid", "LdMRI", "RdMRI")
colnames(dMRIall) = c("eid", "bothdMRI")
colnames(combi_pred) = c("eid", "Lcombi", "Rcombi")
colnames(combiall) = c("eid", "bothcombi")

# merge data frames to match subjects
df01  = (merge(T1_pred, T1all, by = "eid"))
df02 = merge(dMRI_pred, dMRIall, by = "eid")
df03 = merge(combi_pred, combiall, by = "eid")
df04 = merge(df01, df02, by = "eid")
df = merge(df04, df03, by = "eid")

# add assessment centre to data frames
sites = data.frame(demo$eid,demo$site_t3)
colnames(sites) = c("eid","site")
df = merge(df,sites,by = "eid")

## CORRELATE BRAIN AGES AND AGE ####

# show correlation
df %>% dplyr::select(-eid, -site) %>% cor()
res = df %>% dplyr::select(-eid, -site) %>% cor()
breaklist = seq(-1,1, by = 0.1)
my.colors <- (colorRampPalette(colors = c("#E69F00", "darkred"))(length(breaklist)/2))

cors = as.ggplot(pheatmap(res,cluster_rows = F,cluster_cols = F, display_numbers = TRUE, color = my.colors))
# cors + scale_fill_continuous(limits = c(-1,1), breaks = breaklist)

ggsave("/tsd/p33/home/p33-maxk/export/HemiCor.pdf",cors, height = 9, width = 9)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_Plot1.pdf",cors, height = 9, width = 9)


## AGE DIFFERENCES ####


# prep data frames
handedness = handedness %>% dplyr::rename("Handedness"="X1707.0.0") %>% dplyr::select(eid, Handedness)
handedness$Handedness = ordered(handedness$Handedness, levels=c(1,2,3), labels=c("right handed","left handed","ambidextrous"))
df2 = merge(handedness, df, by = "eid")
demo1 = demo %>% dplyr::select(eid, sex, site_t3)
df2 =merge(df2, demo1, by = "eid")

# establish ground truth: age by handedness
df2 %>% group_by(Handedness) %>% summarize(Mean = mean(age), SD = sd(age), N = length(age))
kruskal.test(df2$age ~ df2$Handedness)

#
#
#
# There are age-differences naturally occurring in our sample, so we have to control for age.(Which we would anyways have to do due to the sample-induced age-bias in brain age, but just to mention this here)
#
#
#

############# SOME PREP FOR MODELLING

# look at differences between predictions by hemisphere (for each modality at a time)
## make data frames for each modality
df3.1 = df2 %>% dplyr::select(LT1, RT1) %>% melt()
df3.2 = df2 %>% dplyr::select(LdMRI, RdMRI) %>% melt()
df3.3 = df2 %>% dplyr::select(Lcombi, Rcombi) %>% melt()

## re-add covariates and variable of interest (handedness), hemisphere is indicated in the melted frame as "variable"
df3.1$handedness = c(replicate(2, df2$Handedness))
df3.1$sex = c(replicate(2, df2$sex))
df3.1$age = c(replicate(2, df2$age))
df3.1$site = c(replicate(2, df2$site))
df3.1$eid = c(replicate(2, df2$eid))

df3.2$handedness = c(replicate(2, df2$Handedness))
df3.2$sex = c(replicate(2, df2$sex))
df3.2$age = c(replicate(2, df2$age))
df3.2$site = c(replicate(2, df2$site))
df3.2$eid = c(replicate(2, df2$eid))

df3.3$handedness = c(replicate(2, df2$Handedness))
df3.3$sex = c(replicate(2, df2$sex))
df3.3$age = c(replicate(2, df2$age))
df3.3$site = c(replicate(2, df2$site))
df3.3$eid = c(replicate(2, df2$eid))

# make handedness a factor
df3.1$handedness=as.factor(df3.1$handedness)
df3.2$handedness=as.factor(df3.2$handedness)
df3.3$handedness=as.factor(df3.3$handedness)

# name hemisphere more clearly in the data frames
df3.1 = df3.1 %>% dplyr::rename("hemisphere" = "variable")
df3.2 = df3.2 %>% dplyr::rename("hemisphere" = "variable")
df3.3 = df3.3 %>% dplyr::rename("hemisphere" = "variable")


## SIMPLE BRAIN AGE DIFFERENCES ####
# runs simple linear models to compare hemispheres within each modality while controlling for age
# this is just a preliminary check before employing the mixed model
model01 = lm(value ~ age + hemisphere, data = df3.1)
summary(model01)
model02 = lm(value ~ age + hemisphere, data = df3.2)
summary(model02)
model03 = lm(value ~ age + hemisphere, data = df3.3)
summary(model03)

## HYPOTHESES 1-2: BRAIN AGE DEPENDENCE ON MODALITY AND HEMISPHERE ####
ModHem = rbind(df3.1,df3.2,df3.3)
ModHem$modality = c(replicate(nrow(df3.1),"T1w"),replicate(nrow(df3.2),"dwMRI"),replicate(nrow(df3.3),"multimodal"))
ModHem$hemisphere = c(replicate(nrow(ModHem)/6,"L"),replicate(nrow(ModHem)/6,"R"),replicate(nrow(ModHem)/6,"L"),replicate(nrow(ModHem)/6,"R"),replicate(nrow(ModHem)/6,"L"),replicate(nrow(ModHem)/6,"R"))
modelModHem = lmer(value~hemisphere*modality+sex*age+site+(1|eid), data = (ModHem))
summary(modelModHem)
# sex stratified ####
males = ModHem %>% filter(sex == 1)
females = ModHem %>% filter(sex == 0)
# analyse for females
modelModHem1 = lmer(value~hemisphere*modality+age+(1|site)+(1|eid), data = (females))
summary(modelModHem1)
# analyse for males
modelModHem2 = lmer(value~hemisphere*modality+age+(1|site)+(1|eid), data = (males))
summary(modelModHem2)

## HYPOTHESIS 3: ADDING HANDEDNESS TO THE MODEL ####
# re-estimate models without NA (checking for effect of hemi and interaction effects)
modelModHem = lmer(value~hemisphere*modality+sex*age+(1|site)+(1|eid), data = na.omit(ModHem))
# modelModHemHAND1 = lmer(value~handedness+hemisphere + hemisphere*modality+sex*age+(1|site)+(1|eid), data = na.omit(ModHem))
# modelModHemHAND2 = lmer(value~handedness:hemisphere + hemisphere*modality+sex*age+(1|site)+(1|eid), data = na.omit(ModHem))
modelModHemHAND3 = lmer(value~handedness*hemisphere + hemisphere*modality+sex*age+(1|site)+(1|eid), data = na.omit(ModHem))

# anova(modelModHem,modelModHemHAND1)
# anova(modelModHem,modelModHemHAND2)
anova(modelModHem,modelModHemHAND3)
summary(modelModHemHAND3)

# sex stratification ####
# females
modelModHem1 = lmer(value~hemisphere*modality+age+(1|site)+(1|eid), data = na.omit(females))
modelModHemHAND3_females = lmer(value~handedness*hemisphere + hemisphere*modality+age+(1|site)+(1|eid), data = na.omit(females))
anova(modelModHem1,modelModHemHAND3_females)
summary(modelModHemHAND3_females)
# males
modelModHem2 = lmer(value~hemisphere*modality+age+(1|site)+(1|eid), data = na.omit(males))
modelModHemHAND3_males = lmer(value~handedness*hemisphere + hemisphere*modality+age+(1|site)+(1|eid), data = na.omit(males))
anova(modelModHem2,modelModHemHAND3_males)
summary(modelModHemHAND3_males)

## TESTING ASSOCIATIONS WITH COARIATES ####
# we want diabetics, hypertension, pulse pressure, smoking, vascular diagnosis, waist to hip ratio
vois = vois %>% dplyr::select(-X, -Assessment_centre,-alcohol_drinker)
df[2:11] = scale(df[2:11])
vois$WHR = scale(vois$WHR)
cov_dat = merge(df,vois,by="eid")
# loop over cov_dat frame predicting the different brain ages from the phenotypes of interest and covariates
pheno_vars = colnames(cov_dat[14:ncol(cov_dat)])
pheno_vars = pheno_vars[1:15]
dMRI_vars = colnames(cov_dat[3:11])
RES = list()
RI = list()
for (i in dMRI_vars){
  for (j in pheno_vars){
    f = formula(paste(i,"~",j,"+age*sex+(1|site)"))
    RI[[j]] = lmer(f, data = cov_dat)
  }
  RES[[i]] = RI
}
# get beta values for comparison
Bdf = data.frame(matrix(nrow = length(RES[[1]]),ncol = length(RES)))
Bt = Bdf
Bp = Bdf
B = list()
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bdf[j,i] = summary(RES[[i]][[j]])$coefficients[2,1]
  }
}
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bt[j,i] = summary(RES[[i]][[j]])$coefficients[2,4]
  }
} 
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bp[j,i] = summary(RES[[i]][[j]])$coefficients[2,5]
  }
}
colnames(Bdf) = dMRI_vars
colnames(Bt) = dMRI_vars
colnames(Bp) = dMRI_vars
HealthVars = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
               "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
               "High Cholsterol", "Alcohol", "Smoking")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars

# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_BetaVals.csv")
write.csv(Bt,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_Tvals.csv")
write.csv(Bp,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_Pvals.csv")

Bdfall = Bdf
Btall = Bt
Bpall = Bp

# sex stratification males ####
vois_m = vois %>% dplyr::filter(sex == 1)
cov_dat = merge(df,vois_m,by="eid")
# loop over cov_dat frame predicting the different brain ages from the phenotypes of interest and covariates
pheno_vars = colnames(cov_dat[14:ncol(cov_dat)])
pheno_vars = pheno_vars[1:15]
dMRI_vars = colnames(cov_dat[3:11])
colnames(cov_dat)
RES = list()
RI = list()
for (i in dMRI_vars){
  for (j in pheno_vars){
    f = formula(paste(i,"~",j,"+age+(1|site)"))
    RI[[j]] = lmer(f, data = cov_dat)
  }
  RES[[i]] = RI
}
# get beta values for comparison
Bdf = data.frame(matrix(nrow = length(RES[[1]]),ncol = length(RES)))
Bt = Bdf
Bp = Bdf
B = list()
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bdf[j,i] = summary(RES[[i]][[j]])$coefficients[2,1]
  }
}
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bt[j,i] = summary(RES[[i]][[j]])$coefficients[2,4]
  }
} 
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bp[j,i] = summary(RES[[i]][[j]])$coefficients[2,5]
  }
}
colnames(Bdf) = dMRI_vars
colnames(Bt) = dMRI_vars
colnames(Bp) = dMRI_vars
HealthVars = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
               "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
               "High Cholsterol", "Alcohol", "Smoking")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars
Bdf1 = Bdf
Bt1 = Bt
Bp1 = Bp
# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf1, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_BetaVals.csv")
write.csv(Bt1,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_Tvals.csv")
write.csv(Bp1,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_Pvals.csv")

# sex stratification females ####
vois_f = vois %>% dplyr::filter(sex == 0)
cov_dat = merge(df,vois_f,by="eid")
# loop over cov_dat frame predicting the different brain ages from the phenotypes of interest and covariates
pheno_vars = colnames(cov_dat[14:ncol(cov_dat)])
pheno_vars = pheno_vars[1:15]
dMRI_vars = colnames(cov_dat[3:11])
colnames(cov_dat)
RES = list()
RI = list()
for (i in dMRI_vars){
  for (j in pheno_vars){
    f = formula(paste(i,"~",j,"+age+(1|site)"))
    RI[[j]] = lmer(f, data = cov_dat)
  }
  RES[[i]] = RI
}
# get beta values for comparison
Bdf = data.frame(matrix(nrow = length(RES[[1]]),ncol = length(RES)))
Bt = Bdf
Bp = Bdf
B = list()
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bdf[j,i] = summary(RES[[i]][[j]])$coefficients[2,1]
  }
}
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bt[j,i] = summary(RES[[i]][[j]])$coefficients[2,4]
  }
} 
for (i in 1:length(dMRI_vars)){
  for(j in 1:length(pheno_vars)){
    Bp[j,i] = summary(RES[[i]][[j]])$coefficients[2,5]
  }
}
colnames(Bdf) = dMRI_vars
colnames(Bt) = dMRI_vars
colnames(Bp) = dMRI_vars
HealthVars = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
               "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
               "High Cholsterol", "Alcohol", "Smoking")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars

Bdf0 = Bdf
Bt0 = Bt
Bp0 = Bp

# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf0, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_BetaVals.csv")
write.csv(Bt0,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_Tvals.csv")
write.csv(Bp0,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_Pvals.csv")


# plotting (non sex-stratified) ####
# transform to long
plot_dat = melt(Bdfall)
plot_t = melt(Btall)
plot_p = melt(Bpall)
# add t and p values to the long beta data frame
plot_dat$p = plot_p$value
plot_dat$t = plot_t$value
# adjust p-val
plot_dat$p_adj = p.adjust(plot_dat$p, method = "bonferroni")
pedcolors = c(ifelse(plot_dat$p_adj < .05,"black", "white"))
# round  beta value for visualisation
plot_dat$value = round(plot_dat$value, digits = 2)
plot = ggplot(plot_dat, aes(AssociationVar, variable, fill = value, label = value)) +#col = p, 
  geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(low = "white", mid = "azure1",high = "lightblue") +
  #scale_color_gradient2(low = "red", mid = "orange", high = "white") +
  ylab("Brain Age")  + labs(fill = "Beta Value") + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  xlab("Phenotype") #+ geom_tile(data = frames, fill = "black")
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot.pdf",plot, height =8, width = 10)


# plotting (males only) ####
# transform to long
plot_dat = melt(Bdf1)
plot_t = melt(Bt1)
plot_p = melt(Bp1)
# add t and p values to the long beta data frame
plot_dat$p = plot_p$value
plot_dat$t = plot_t$value
# adjust p-val
plot_dat$p_adj = p.adjust(plot_dat$p, method = "bonferroni")
pedcolors = c(ifelse(plot_dat$p_adj < .05,"black", "white"))
# round  beta value for visualisation
plot_dat$value = round(plot_dat$value, digits = 2)
plot1 = ggplot(plot_dat, aes(AssociationVar, variable, fill = value, label = value)) +#col = p, 
  geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(low = "white", mid = "azure1",high = "lightblue") +
  #scale_color_gradient2(low = "red", mid = "orange", high = "white") +
  ylab("Brain Age")  + labs(fill = "Beta Value") + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  xlab("Phenotype") #+ geom_tile(data = frames, fill = "black")
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_males.pdf",plot, height =8, width = 10)


# plotting (females only) ####
# transform to long
plot_dat = melt(Bdf0)
plot_t = melt(Bt0)
plot_p = melt(Bp0)
# add t and p values to the long beta data frame
plot_dat$p = plot_p$value
plot_dat$t = plot_t$value
# adjust p-val
plot_dat$p_adj = p.adjust(plot_dat$p, method = "bonferroni")
pedcolors = c(ifelse(plot_dat$p_adj < .05,"black", "white"))
# round  beta value for visualisation
plot_dat$value = round(plot_dat$value, digits = 2)
plot2 = ggplot(plot_dat, aes(AssociationVar, variable, fill = value, label = value)) +#col = p, 
  geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(low = "white",high = "lightblue") +
  #scale_color_gradient2(low = "red", mid = "orange", high = "white") +
  ylab("Brain Age")  + labs(fill = "Beta Value") + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  xlab("Phenotype") #+ geom_tile(data = frames, fill = "black")
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_females.pdf",plot, height =8, width = 10)

# save both sex plots together
sex_plot = ggarrange(plot1, plot2, labels = c("a", "b"), common.legend = T)

#ggsave("/tsd/p33/home/p33-maxk/export/sex_plot.pdf", sex_plot, height =8, width = 17)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_sex_stratification.pdf",sex_plot, height =8, width = 17)
