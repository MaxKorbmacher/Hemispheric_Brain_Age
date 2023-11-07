# test laterality index of brain age
# Max Korbmacher, 14.08.2023
# 

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
vois = read.csv("/cluster/projects/p33/users/maxk/UKB/BP_BAG/data/demo_clean.csv")


# put predictions from different hemispheres into one data set
T1_pred = data.frame(LT1$eid,LT1$age, LT1$pred_age_LT1, RT1$pred_age_RT1)
dMRI_pred = data.frame(LdMRI$eid, LdMRI$pred_age_LdMRI,RdMRI$pred_age_RdMRI)
combi_pred = data.frame(Lcombi$eid, Lcombi$pred_age_Lcombi,Rcombi$pred_age_Rcombi)

# change column names
colnames(T1_pred) = c("eid","age", "LT1", "RT1")
colnames(dMRI_pred) = c("eid", "LdMRI", "RdMRI")
colnames(combi_pred) = c("eid", "Lcombi", "Rcombi")

# prepare data frame to loop over for the estimation of corrected brain ages
df01  = (merge(T1_pred, dMRI_pred, by = "eid"))
df02  = (merge(df01, combi_pred, by = "eid"))
df02 = df02 %>% dplyr::relocate(c(eid,age),.after=last_col())

# estimate LI for uncorrected brain ages
T1_pred$T1LI = abs((T1_pred$LT1 - T1_pred$RT1)/(T1_pred$LT1 + T1_pred$RT1))
dMRI_pred$dMRILI = abs((dMRI_pred$LdMRI - dMRI_pred$RdMRI)/(dMRI_pred$LdMRI + dMRI_pred$RdMRI))
combi_pred$combiLI = abs((combi_pred$Lcombi - combi_pred$Rcombi)/(combi_pred$Lcombi + combi_pred$Rcombi))

# merge data frames of uncorrected brain ages to match subjects
df01  = (merge(T1_pred, dMRI_pred, by = "eid"))
df = merge(combi_pred, df01, by = "eid")

# add assessment centre and sex to data frames
sites = data.frame(demo$eid,demo$site_t3, demo$sex)
colnames(sites) = c("eid","site", "sex")
df = merge(df,sites,by = "eid")

# select only variables of interest
df = df %>% dplyr::select(combiLI, T1LI,dMRILI,age,eid, sex, site)

# standardize values
df[1:4] = df[1:4]%>%mutate_if(is.numeric,scale)

# estimate associations between brain age LIs and age######
brainages = c("combiLI", "T1LI","dMRILI")
for (i in brainages){
  # run mixed effects model
  f = formula(paste("age ~", i, "+ sex + (1|site)"))
  model1 = lmer(f,data = df)
  # run control model
  f2 = formula(paste("age ~ sex + (1|site)"))
  model2 = lmer(f2,data = df)
  # print stats for the relationship of LI brain age and age
  print(paste("Adjusted association for ",i))
  print(summary(model1)$coefficients[2,])
  # print crude correlation and p-value
  f3 = formula(paste("age ~",i))
  model3 = lm(f3,data = df)
  print(paste("Crude correlation for ",i))
  print(model3$coefficients[2])
  print(paste("P-value for crude association for ",i))
  print(summary(model3)$coefficients[2,4])
  # print likelihood ratio test outcome
  print(paste("LR Test Outcome for",i))
  print(anova(model2,model1))
  # mark end of one brain age estimation
  print("###########################")
}
# repeating the analyses for males only ####
df1 = df %>% filter(sex == 1)
for (i in brainages){
  # run mixed effects model
  f = formula(paste("age ~", i, "+ (1|site)"))
  model1 = lmer(f,data = df1)
  # run control model
  f2 = formula(paste("age ~ (1|site)"))
  model2 = lmer(f2,data = df1)
  # print stats for the relationship of LI brain age and age
  print(paste("Adjusted association for ",i))
  print(summary(model1)$coefficients[2,])
  # print crude correlation and p-value
  f3 = formula(paste("age ~",i))
  model3 = lm(f3,data = df1)
  print(paste("Crude correlation for ",i))
  print(model3$coefficients[2])
  print(paste("P-value for crude association for ",i))
  print(summary(model3)$coefficients[2,4])
  # print likelihood ratio test outcome
  print(paste("LR Test Outcome for",i))
  print(anova(model2,model1))
  # mark end of one brain age estimation
  print("###########################")
}
# repeating the analyses for females only ####
df0 = df %>% filter(sex == 0)
for (i in brainages){
  # run mixed effects model
  f = formula(paste("age ~", i, "+ (1|site)"))
  model1 = lmer(f,data = df0)
  # run control model
  f2 = formula(paste("age ~ (1|site)"))
  model2 = lmer(f2,data = df0)
  # print stats for the relationship of LI brain age and age
  print(paste("Adjusted association for ",i))
  print(summary(model1)$coefficients[2,])
  # print crude correlation and p-value
  f3 = formula(paste("age ~",i))
  model3 = lm(f3,data = df0)
  print(paste("Crude correlation for ",i))
  print(model3$coefficients[2])
  print(paste("P-value for crude association for ",i))
  print(summary(model3)$coefficients[2,4])
  # print likelihood ratio test outcome
  print(paste("LR Test Outcome for",i))
  print(anova(model2,model1))
  # mark end of one brain age estimation
  print("###########################")
}
