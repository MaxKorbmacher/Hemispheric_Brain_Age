# exploring hemispheric age predictions
# updated: 02.11.2023 adding sex stratification

## PREP ####

# packages
#install.packages("afex")
#install.packages("tidyquant")
library(tidyquant)
#install.packages("FSA")
library(FSA)
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
library(gridExtra)
#install.packages("marginaleffects")
library(marginaleffects)
library(lmerTest)

# read in data (both sexes)
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
# males
LT1_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LT1_males.csv")
RT1_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RT1_males.csv")
LdMRI_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LdMRI_males.csv")
RdMRI_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RdMRI_males.csv")
Lcombi_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_Lcombi_males.csv")
Rcombi_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_Rcombi_males.csv")
bothT1_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_T1_males.csv")
bothdMRI_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_dMRI_males.csv")
bothcombi_males = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_combi_males.csv")
# females
LT1_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LT1_females.csv")
RT1_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RT1_females.csv")
LdMRI_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_LdMRI_females.csv")
RdMRI_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_RdMRI_females.csv")
Lcombi_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/L/Brainage_Lcombi_females.csv")
Rcombi_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/R/Brainage_Rcombi_females.csv")
bothT1_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_T1_females.csv")
bothdMRI_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_dMRI_females.csv")
bothcombi_females = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/Brainage_combi_females.csv")
# variables of interest
handedness = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/analyses/handedness.csv")
vois = read.csv("/cluster/projects/p33/users/maxk/UKB/BP_BAG/data/demo_clean.csv")

# put predictions from different hemispheres into one data set
## both sexes
T1_pred = data.frame(LT1$eid,LT1$age, LT1$pred_age_LT1, RT1$pred_age_RT1)
T1all = data.frame(bothT1$eid,bothT1$pred_age_T1)
dMRI_pred = data.frame(LdMRI$eid, LdMRI$pred_age_LdMRI,RdMRI$pred_age_RdMRI)
dMRIall = data.frame(bothdMRI$eid,bothdMRI$pred_age_dwMRI)
combi_pred = data.frame(Lcombi$eid, Lcombi$pred_age_Lcombi,Rcombi$pred_age_Rcombi)
combiall = data.frame(bothcombi$eid,bothcombi$pred_age_combi)
## males
T1_pred_males = data.frame(LT1_males$eid,LT1_males$age, LT1_males$pred_age_LT1, RT1_males$pred_age_RT1)
T1all_males = data.frame(bothT1_males$eid,bothT1_males$pred_age_T1)
dMRI_pred_males = data.frame(LdMRI_males$eid, LdMRI_males$pred_age_LdMRI,RdMRI_males$pred_age_RdMRI)
dMRIall_males = data.frame(bothdMRI_males$eid,bothdMRI_males$pred_age_dMRI)
combi_pred_males = data.frame(Lcombi_males$eid, Lcombi_males$pred_age_Lcombi,Rcombi_males$pred_age_Rcombi)
combiall_males = data.frame(bothcombi_males$eid,bothcombi_males$pred_age_combi)
## females
T1_pred_females = data.frame(LT1_females$eid,LT1_females$age, LT1_females$pred_age_LT1, RT1_females$pred_age_RT1)
T1all_females = data.frame(bothT1_females$eid,bothT1_females$pred_age_T1)
dMRI_pred_females = data.frame(LdMRI_females$eid, LdMRI_females$pred_age_LdMRI,RdMRI_females$pred_age_RdMRI)
dMRIall_females = data.frame(bothdMRI_females$eid,bothdMRI_females$pred_age_dMRI)
combi_pred_females = data.frame(Lcombi_females$eid, Lcombi_females$pred_age_Lcombi,Rcombi_females$pred_age_Rcombi)
combiall_females = data.frame(bothcombi_females$eid,bothcombi_females$pred_age_combi)

# change column names (both sexes)
colnames(T1_pred) = c("eid","age", "LT1", "RT1")
colnames(T1all) = c("eid", "bothT1")
colnames(dMRI_pred) = c("eid", "LdMRI", "RdMRI")
colnames(dMRIall) = c("eid", "bothdMRI")
colnames(combi_pred) = c("eid", "Lcombi", "Rcombi")
colnames(combiall) = c("eid", "bothcombi")
# change column names (males)
colnames(T1_pred_males) = colnames(T1_pred)
colnames(T1all_males) = colnames(T1all) 
colnames(dMRI_pred_males) = colnames(dMRI_pred) 
colnames(dMRIall_males) = colnames(dMRIall)
colnames(combi_pred_males) = colnames(combi_pred)
colnames(combiall_males) = colnames(combiall)
# change column names (females)
colnames(T1_pred_females) = colnames(T1_pred)
colnames(T1all_females) = colnames(T1all) 
colnames(dMRI_pred_females) = colnames(dMRI_pred) 
colnames(dMRIall_females) = colnames(dMRIall)
colnames(combi_pred_females) = colnames(combi_pred)
colnames(combiall_females) = colnames(combiall)

# merge data frames to match subjects (both sexes)
df01  = (merge(T1_pred, T1all, by = "eid"))
df02 = merge(dMRI_pred, dMRIall, by = "eid")
df03 = merge(combi_pred, combiall, by = "eid")
df04 = merge(df01, df02, by = "eid")
df = merge(df04, df03, by = "eid")

# merge data frames to match subjects (males)
df01_males  = (merge(T1_pred_males, T1all_males, by = "eid"))
df02_males = merge(dMRI_pred_males, dMRIall_males, by = "eid")
df03_males = merge(combi_pred_males, combiall_males, by = "eid")
df04_males = merge(df01_males, df02_males, by = "eid")
df_males = merge(df04_males, df03_males, by = "eid")

# merge data frames to match subjects (females)
df01_males  = (merge(T1_pred_females, T1all_females, by = "eid"))
df02_males = merge(dMRI_pred_females, dMRIall_females, by = "eid")
df03_males = merge(combi_pred_females, combiall_females, by = "eid")
df04_males = merge(df01_males, df02_males, by = "eid")
df_females = merge(df04_males, df03_males, by = "eid")

# add assessment centre to data frames
sites = data.frame(demo$eid,demo$site_t3)
colnames(sites) = c("eid","site")
df = merge(df,sites,by = "eid")
df_males = merge(df_males,sites,by = "eid")
df_females = merge(df_females,sites,by = "eid")

## CORRELATE BRAIN AGES AND AGE ####

# show correlation
df %>% dplyr::select(-eid, -site) %>% cor()
res = df %>% dplyr::select(-eid, -site)
names(res) = c("Age", "T1: L", "T1: R", "T1: LR", "diffusion: L", "diffusion: R", "diffusion: LR", "multimodal: L", "multimodal: R", "multimodal: LR")
res = res %>% cor()
breaklist = seq(-1,1, by = 0.1)
my.colors <- (colorRampPalette(colors = c("#E69F00", "darkred"))(length(breaklist)/2))

cors = as.ggplot(pheatmap(res,cluster_rows = F,cluster_cols = F, display_numbers = TRUE, color = my.colors))
# cors + scale_fill_continuous(limits = c(-1,1), breaks = breaklist)

ggsave("/tsd/p33/home/p33-maxk/export/HemiCor.pdf",cors, height = 9, width = 9)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_Plot1.pdf",cors, height = 9, width = 9)

# do the same for sex specific models (trained only on males' or females' data)
res1 = df_males %>% dplyr::select(-eid, -site)
res2 = df_females %>% dplyr::select(-eid, -site)
names(res1) = c("Age", "T1: L", "T1: R", "T1: LR", "diffusion: L", "diffusion: R", "diffusion: LR", "multimodal: L", "multimodal: R", "multimodal: LR")
names(res2) = c("Age", "T1: L", "T1: R", "T1: LR", "diffusion: L", "diffusion: R", "diffusion: LR", "multimodal: L", "multimodal: R", "multimodal: LR")
res1 = res1 %>% cor()
res2 = res2 %>% cor()
p1=pheatmap(res1,cluster_rows = F,cluster_cols = F, display_numbers = TRUE, color = my.colors, main = "Males")
p2=pheatmap(res2,cluster_rows = F,cluster_cols = F, display_numbers = TRUE, color = my.colors, main = "Females")
plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]
a = grid.arrange(grobs=plot_list, ncol=2)
ggsave("/tsd/p33/home/p33-maxk/export/HemiCor_sexes.pdf",a, height = 6, width = 14)
ggsave("/cluster/projects/p33/users/maxk/UKB/export//Hemi_NEW_sex_stratefied_brain_age_correlations.pdf",a, height = 6, width = 14)


## AGE DIFFERENCES ####


# prep data frames
handedness = handedness %>% dplyr::rename("Handedness"="X1707.0.0") %>% dplyr::select(eid, Handedness)
handedness$Handedness = ordered(handedness$Handedness, levels=c(1,2,3), labels=c("right handed","left handed","ambidextrous"))
df2 = merge(handedness, df, by = "eid")
df2_males = merge(handedness, df_males, by = "eid")
df2_females = merge(handedness, df_females, by = "eid")
demo1 = demo %>% dplyr::select(eid, sex, site_t3)
df2 =merge(df2, demo1, by = "eid")
df2_males =merge(df2_males, demo1, by = "eid")
df2_females =merge(df2_females, demo1, by = "eid")

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
# males 
df3.1_males = df2_males %>% dplyr::select(LT1, RT1) %>% melt()
df3.2_males = df2_males %>% dplyr::select(LdMRI, RdMRI) %>% melt()
df3.3_males = df2_males %>% dplyr::select(Lcombi, Rcombi) %>% melt()
# females 
df3.1_females = df2_females %>% dplyr::select(LT1, RT1) %>% melt()
df3.2_females = df2_females %>% dplyr::select(LdMRI, RdMRI) %>% melt()
df3.3_females = df2_females %>% dplyr::select(Lcombi, Rcombi) %>% melt()


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

# do the same for male and female specific brain ages
df3.1_males$handedness = c(replicate(2, df2_males$Handedness))
df3.1_males$sex = c(replicate(2, df2_males$sex))
df3.1_males$age = c(replicate(2, df2_males$age))
df3.1_males$site = c(replicate(2, df2_males$site))
df3.1_males$eid = c(replicate(2, df2_males$eid))

df3.2_males$handedness = c(replicate(2, df2_males$Handedness))
df3.2_males$sex = c(replicate(2, df2_males$sex))
df3.2_males$age = c(replicate(2, df2_males$age))
df3.2_males$site = c(replicate(2, df2_males$site))
df3.2_males$eid = c(replicate(2, df2_males$eid))

df3.3_males$handedness = c(replicate(2, df2_males$Handedness))
df3.3_males$sex = c(replicate(2, df2_males$sex))
df3.3_males$age = c(replicate(2, df2_males$age))
df3.3_males$site = c(replicate(2, df2_males$site))
df3.3_males$eid = c(replicate(2, df2_males$eid))

# make handedness a factor
df3.1_males$handedness=as.factor(df3.1_males$handedness)
df3.2_males$handedness=as.factor(df3.2_males$handedness)
df3.3_males$handedness=as.factor(df3.3_males$handedness)

# name hemisphere more clearly in the data frames
df3.1_males = df3.1_males %>% dplyr::rename("hemisphere" = "variable")
df3.2_males = df3.2_males %>% dplyr::rename("hemisphere" = "variable")
df3.3_males = df3.3_males %>% dplyr::rename("hemisphere" = "variable")

# ... for females
df3.1_females$handedness = c(replicate(2, df2_females$Handedness))
df3.1_females$sex = c(replicate(2, df2_females$sex))
df3.1_females$age = c(replicate(2, df2_females$age))
df3.1_females$site = c(replicate(2, df2_females$site))
df3.1_females$eid = c(replicate(2, df2_females$eid))

df3.2_females$handedness = c(replicate(2, df2_females$Handedness))
df3.2_females$sex = c(replicate(2, df2_females$sex))
df3.2_females$age = c(replicate(2, df2_females$age))
df3.2_females$site = c(replicate(2, df2_females$site))
df3.2_females$eid = c(replicate(2, df2_females$eid))

df3.3_females$handedness = c(replicate(2, df2_females$Handedness))
df3.3_females$sex = c(replicate(2, df2_females$sex))
df3.3_females$age = c(replicate(2, df2_females$age))
df3.3_females$site = c(replicate(2, df2_females$site))
df3.3_females$eid = c(replicate(2, df2_females$eid))

# make handedness a factor
df3.1_females$handedness=as.factor(df3.1_females$handedness)
df3.2_females$handedness=as.factor(df3.2_females$handedness)
df3.3_females$handedness=as.factor(df3.3_females$handedness)

# name hemisphere more clearly in the data frames
df3.1_females = df3.1_females %>% dplyr::rename("hemisphere" = "variable")
df3.2_females = df3.2_females %>% dplyr::rename("hemisphere" = "variable")
df3.3_females = df3.3_females %>% dplyr::rename("hemisphere" = "variable")

## SIMPLE BRAIN AGE DIFFERENCES ####
# runs simple linear models to compare hemispheres within each modality while controlling for age
# this is just a preliminary check before employing the mixed model
model01 = lm(value ~ age + hemisphere, data = df3.1) # T1
summary(model01)
model02 = lm(value ~ age + hemisphere, data = df3.2) #dMRI
summary(model02)
model03 = lm(value ~ age + hemisphere, data = df3.3) # multimodal
summary(model03)

# estimate differences when considering females or males brain ages separately or together, but controlling for age and scanner site
# T1w
merger = df3.1 %>% filter(sex == 1)
malesT1 = rbind(merger, df3.1_males)
malesT1$model = c(replicate(nrow(merger), "MalesBoth"), replicate(nrow(merger), "MalesMales"))
merger = df3.1 %>% filter(sex == 0)
femalesT1 = rbind(merger, df3.1_females)
femalesT1$model = c(replicate(nrow(merger), "FemalesBoth"), replicate(nrow(merger), "FemalesFemales"))
T1_final = rbind(malesT1, femalesT1)
# dMRI
merger = df3.2 %>% filter(sex == 1)
malesdMRI = rbind(merger, df3.2_males)
malesdMRI$model = c(replicate(nrow(merger), "MalesBoth"), replicate(nrow(merger), "MalesMales"))
merger = df3.2 %>% filter(sex == 0)
femalesdMRI = rbind(merger, df3.2_females)
femalesdMRI$model = c(replicate(nrow(merger), "FemalesBoth"), replicate(nrow(merger), "FemalesFemales"))
dMRI_final = rbind(malesdMRI, femalesdMRI)
# multimodal
merger = df3.3 %>% filter(sex == 1)
malesmulti = rbind(merger, df3.3_males)
malesmulti$model = c(replicate(nrow(merger), "MalesBoth"), replicate(nrow(merger), "MalesMales"))
merger = df3.3 %>% filter(sex == 0)
femalesmulti = rbind(merger, df3.3_females)
femalesmulti$model = c(replicate(nrow(merger), "FemalesBoth"), replicate(nrow(merger), "FemalesFemales"))
multi_final = rbind(malesmulti, femalesmulti)

# now check whether estimates differ
final_data = rbind(T1_final,dMRI_final,multi_final)
final_data$model = factor(final_data$model)
final_data$sex = factor(final_data$sex)
final_data$modality = c(replicate(nrow(T1_final),"T1w"),replicate(nrow(dMRI_final),"dwMRI"),replicate(nrow(multi_final),"multimodal"))

# lmer
all_mod = lmer(value~hemisphere*modality*sex+age+model+(1|site)+(1|eid), data = final_data)

summary(all_mod)

# check whether marginal estimates differ between models (sex specific vs both sexes)
## for females
avg_comparisons(all_mod,variables = list(model = c("FemalesBoth", "FemalesFemales")))
## for males
avg_comparisons(all_mod,variables = list(model = c("MalesBoth", "MalesMales")))

# sex effect
avg_comparisons(all_mod,variables = list(sex = c("0","1")))

# hemisphere effect
avg_comparisons(all_mod,variables = list(hemisphere = c("LT1","RT1")))
avg_comparisons(all_mod,variables = list(hemisphere = c("LdMRI","RdMRI")))
avg_comparisons(all_mod,variables = list(hemisphere = c("Lcombi","Rcombi")))

# handedness effect
avg_comparisons(all_mod,variables = list(handedness = c("ambidextrous","left handed", "right handed")))

final_data$prediction = predict(all_mod)
kruskal.test(prediction ~ handedness, data = final_data)
final_data %>% group_by(handedness) %>% summarize(Mean = mean(prediction), SD = sd(prediction), N = length(prediction))
# post hoc test showing handedness differences
dunnTest(value ~ handedness,
              data=final_data,
              method="bh")

# show interaction effects
avg_comparisons(all_mod,variables = "sex", by = "hemisphere")

avg_comparisons(all_mod,variables = list(hemisphere = c("LT1","RT1")), by = "sex")
avg_comparisons(all_mod,variables = list(hemisphere = c("LdMRI","RdMRI")), by = "sex")
avg_comparisons(all_mod,variables = list(hemisphere = c("Lcombi","Rcombi")), by = "sex")

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
# follow-up investigating sex-dependent hemispheric effects
modelModHemHAND4 = lmer(value~handedness*hemisphere + hemisphere*modality+hemisphere:sex+modality:sex+sex*age+(1|site)+(1|eid), data = na.omit(ModHem))
summary(modelModHemHAND4)

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


# PLOTS WHICH ALSO CONSIDER SEX-SPECIFIC TRAINING ####
T1_pred_males = data.frame(LT1_males$eid,LT1_males$age, LT1_males$pred_age_LT1, RT1_males$pred_age_RT1, bothT1_males$pred_age_T1_males)
T1_pred_females = data.frame(LT1_females$eid,LT1_females$age, LT1_females$pred_age_LT1, RT1_females$pred_age_RT1, bothT1_females$pred_age_T1_females)
colnames(T1_pred_females) = c("eid","age", "LT1", "RT1", "bothT1")
colnames(T1_pred_males) = colnames(T1_pred_females)
T1_pred = rbind(T1_pred_males, T1_pred_females)
T1_pred$sex = c(replicate(nrow(T1_pred_males), "male"), replicate(nrow(T1_pred_females), "female"))
# dMRI
dMRI_pred_males = data.frame(LdMRI_males$eid,LdMRI_males$age, LdMRI_males$pred_age_LdMRI, RdMRI_males$pred_age_RdMRI, bothdMRI_males$pred_age_dMRI_males)
dMRI_pred_females = data.frame(LdMRI_females$eid,LdMRI_females$age, LdMRI_females$pred_age_LdMRI, RdMRI_females$pred_age_RdMRI, bothdMRI_females$pred_age_dMRI_females)
colnames(dMRI_pred_females) = c("eid","age", "LdMRI", "RdMRI", "bothdMRI")
colnames(dMRI_pred_males) = c("eid","age", "LdMRI", "RdMRI", "bothdMRI")
dMRI_pred = rbind(dMRI_pred_males, dMRI_pred_females)
dMRI_pred$sex = c(replicate(nrow(dMRI_pred_males), "male"), replicate(nrow(dMRI_pred_females), "female"))
# multimodal
combi_pred_males = data.frame(Lcombi_males$eid,Lcombi_males$age, Lcombi_males$pred_age_Lcombi, Rcombi_males$pred_age_Rcombi, bothcombi_males$pred_age_combi_males)
combi_pred_females = data.frame(Lcombi_females$eid,Lcombi_females$age, Lcombi_females$pred_age_Lcombi, Rcombi_females$pred_age_Rcombi, bothcombi_females$pred_age_combi_females)
colnames(combi_pred_females) = c("eid","age", "Lcombi", "Rcombi", "bothcombi")
colnames(combi_pred_males) = c("eid","age", "Lcombi", "Rcombi", "bothcombi")
combi_pred = rbind(combi_pred_males, combi_pred_females)
combi_pred$sex = c(replicate(nrow(combi_pred_males), "male"), replicate(nrow(combi_pred_females), "female"))


# prepare data frame to loop over for the estimation of corrected brain ages
df01  = (merge(T1_pred, dMRI_pred, by = "eid"))
df02  = (merge(df01, combi_pred, by = "eid"))
df02 = df02 %>% dplyr::relocate(c(age,eid,sex),.after=last_col())
df = df02 %>% select(-age.x, -age.y, -sex.y, -sex.x)
df[1:9] = scale(df[1:9])
vois$WHR = scale(vois$WHR)
cov_dat = merge(df,vois,by="eid")
cov_dat$sex = cov_dat$sex.x
cov_dat = cov_dat%>% select(-sex.y, -sex.x)
sites = data.frame(demo$eid,demo$site_t3)
colnames(sites) = c("eid","site")
cov_dat = merge(cov_dat,sites,by = "eid")

# loop over cov_dat frame predicting the different brain ages from the phenotypes of interest and covariates
pheno_vars = colnames(cov_dat[12:27])
dMRI_vars = colnames(cov_dat[2:10])
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
               "High Cholsterol", "Alcohol", "Smoking", "Ethnicity")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars

# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_BetaVals_sex_specific_models.csv")
write.csv(Bt,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_Tvals_sex_specific_models.csv")
write.csv(Bp,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_Pvals_sex_specific_models.csv")

Bdfall = Bdf
Btall = Bt
Bpall = Bp

# sex stratification males ####
vois_m = vois %>% dplyr::filter(sex == 1)
cov_dat = merge(df,vois_m,by="eid")
cov_dat = merge(cov_dat,sites,by = "eid")
# loop over cov_dat frame predicting the different brain ages from the phenotypes of interest and covariates
# pheno_vars = colnames(cov_dat[14:ncol(cov_dat)])
# pheno_vars = pheno_vars[1:15]
# dMRI_vars = colnames(cov_dat[3:11])
# colnames(cov_dat)
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
# HealthVars = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
#                "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
#                "High Cholsterol", "Alcohol", "Smoking")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars
Bdf1 = Bdf
Bt1 = Bt
Bp1 = Bp
# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf1, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_BetaVals_sex_specific_models.csv")
write.csv(Bt1,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_Tvals_sex_specific_models.csv")
write.csv(Bp1,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_male_pheno_Pvals_sex_specific_models.csv")

# sex stratification females ####
vois_f = vois %>% dplyr::filter(sex == 0)
cov_dat = merge(df,vois_f,by="eid")
cov_dat = merge(cov_dat,sites,by = "eid")

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
# HealthVars = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
#                "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
#                "High Cholsterol", "Alcohol", "Smoking")
Bdf$AssociationVar = HealthVars
Bt$AssociationVar = HealthVars
Bp$AssociationVar = HealthVars

Bdf0 = Bdf
Bt0 = Bt
Bp0 = Bp

# write results as csv files (t-vals, p-vals, beta values)
write.csv(Bdf0, "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_BetaVals_sex_specific_models.csv")
write.csv(Bt0,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_Tvals_sex_specific_models.csv")
write.csv(Bp0,"/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_female_pheno_Pvals_sex_specific_models.csv")


# plotting (non sex-stratified) ####
# transform to long
Bdfall = Bdfall[1:15,]
Btall = Btall[1:15,]
Bpall = Bpall[1:15,]
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
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_sex_specific_models.pdf",plot, height =8, width = 10)
ggsave("/tsd/p33/home/p33-maxk/export/Hemi_NEW_pheno_plot_sex_specific_models.pdf", plot, height =8, width = 10)


# plotting (males only) ####
# transform to long
Bdf1 = Bdf1[1:15,]
Bt1 = Bt1[1:15,]
Bp1 = Bp1[1:15,]
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
Asso = plot_dat$AssociationVar
plot_dat$value = round(plot_dat$value, digits = 2)
plot1 = ggplot(plot_dat, aes(AssociationVar, variable, fill = value, label = value)) +#col = p, 
  geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
  geom_text(col = "black") +
  theme_minimal() +
  scale_fill_gradient2(low = "white", mid = "azure1",high = "lightblue") +
  #scale_color_gradient2(low = "red", mid = "orange", high = "white") +
  ylab("Brain Age")  + labs(fill = "Beta Value") + theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  xlab("Phenotype") #+ geom_tile(data = frames, fill = "black")
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_males_sex_specific_models.pdf",plot1, height =8, width = 10)
ggsave("/tsd/p33/home/p33-maxk/export/Hemi_NEW_pheno_plot_males_sex_specific_models.pdf", plot1, height =8, width = 10)


# plotting (females only) ####
# transform to long
Bdf0$AssociationVar = c("Waist Circumference", "Hip Circumference", "Height","Smoking Status", "BMI", "Weight", "WHR", 
                                     "Diastolic Blood Pressure", "Systolic Blood Pressure", "Pulse Pressure", "Diabetis", "Hypertension",
                                     "High Cholsterol", "Alcohol", "Smoking")
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
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_females_sex_specific_models.pdf",plot2, height =8, width = 10)

# save both sex plots together
sex_plot = ggarrange(plot1, plot2, labels = c("a", "b"), common.legend = T)

#ggsave("/tsd/p33/home/p33-maxk/export/sex_plot.pdf", sex_plot, height =8, width = 17)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_pheno_plot_sex_stratification_sex_specific_models.pdf",sex_plot, height =8, width = 17)
ggsave("/tsd/p33/home/p33-maxk/export/Hemi_NEW_pheno_plot_sex_stratification_sex_specific_models.pdf",sex_plot, height =8, width = 17)
