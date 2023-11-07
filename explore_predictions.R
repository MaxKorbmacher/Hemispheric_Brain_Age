# exploring hemispheric age predictions


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

handedness = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/analyses/handedness.csv")

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

## CORRELATE ####

# show correlation
df %>% dplyr::select(-eid, -site) %>% cor()
res = df %>% dplyr::select(-eid, -site) %>% cor()
cors = as.ggplot(pheatmap(res,cluster_rows = F,cluster_cols = F, display_numbers = TRUE))
#ggsave("/tsd/p33/home/p33-maxk/export/HemiCor.pdf",cors, height = 9, width = 9)
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
# There are age-differences naturally occurring in our sample, so we have to control for age.(Which we would anyways have to do when treating brain age, but just saying)
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

## HYPOTHESIS 1: BRAIN AGE DEPENDENCE ON MODALITY AND HEMISPHERE ####
ModHem = rbind(df3.1,df3.2)
ModHem$modality = c(replicate(nrow(df3.1),"T1w"),replicate(nrow(df3.2),"dwMRI"))
ModHem$hemisphere = c(replicate(nrow(ModHem)/4,"L"),replicate(nrow(ModHem)/4,"R"),replicate(nrow(ModHem)/4,"L"),replicate(nrow(ModHem)/4,"R"))
modelModHem = lmer(value~hemisphere*modality+sex*age+site+(1|eid), data = (ModHem))
summary(modelModHem)


## HANDEDNESS & HEMISPHERIC DIFFERENCES (MIXED MODELS) ####

# run the three models comparing hemispheres and handedness
T1_null_model = lmer(value~sex*age+site+(1|eid), data = na.omit(df3.1))
T1_model = lmer(value~hemisphere*handedness+sex*age+site+(1|eid), data = na.omit(df3.1))
dMRI_null_model = lmer(value~sex*age+site+(1|eid), data = na.omit(df3.2))
dMRI_model = lmer(value~hemisphere*handedness+sex*age+site+(1|eid), data = na.omit(df3.2))
combi_null_model = lmer(value~sex*age+site+(1|eid), data = na.omit(df3.3))
combi_model = lmer(value~hemisphere*handedness+sex*age+site+(1|eid), data = na.omit(df3.3))

# model differences?
anova(T1_null_model, T1_model)
anova(dMRI_null_model, dMRI_model)
anova(combi_null_model, combi_model)

# model summaries, if wanting to look at details
#summary(T1_model)
#summary(dMRI_model)
#summary(combi_model)
# can also extract a set of parameters, if wanted
#c(r.squaredGLMM(T1_model),AIC(T1_model),(sigma(T1_model))^2,summary(T1_model)$logLik[1])
# for example check the difference in R2 from null to full model from the multimodal brain age
r.squaredGLMM(combi_model)
r.squaredGLMM(combi_null_model)
r.squaredGLMM(combi_model) - r.squaredGLMM(combi_null_model)
# prep beta plot
T1_dat = data.frame(summary(T1_model)$coefficients)
T1_dat$Variables = rownames(T1_dat)
T1_dat
T1_dat = T1_dat[-c(1,5,6,9),] # delete unwanted rows containing beta values for intercept, sex*age (incl interaction term)
T1_dat$Variables = c("R over L Hemisphere", "L handedness over ambidexterity", "R handedness over ambidexterity", "R-Hemisphere-L-Handedness Interaction", "R-Hemisphere-R-Handedness Interaction")
dMRI_dat = data.frame(summary(dMRI_model)$coefficients)
dMRI_dat$Variables = rownames(dMRI_dat)
dMRI_dat = dMRI_dat[-c(1,5,6,9),] # delete unwanted rows containing beta values for intercept, sex*age (incl interaction term)
dMRI_dat$Variables = c("R over L Hemisphere", "L handedness over ambidexterity", "R handedness over ambidexterity", "R-Hemisphere-L-Handedness Interaction", "R-Hemisphere-R-Handedness Interaction")
combi_dat = data.frame(summary(combi_model)$coefficients)
combi_dat$Variables = rownames(combi_dat)
combi_dat = combi_dat[-c(1,5,6,9),] # delete unwanted rows containing beta values for intercept, sex*age (incl interaction term)
combi_dat$Variables = c("R over L Hemisphere", "L handedness over ambidexterity", "R handedness over ambidexterity", "R-Hemisphere-L-Handedness Interaction", "R-Hemisphere-R-Handedness Interaction")
plot_dat = rbind(T1_dat,dMRI_dat,combi_dat)
plot_dat$Metric = c(replicate(nrow(T1_dat), "T1w"),replicate(nrow(dMRI_dat), "dMRI"),replicate(nrow(combi_dat), "T1&dMRI"))
plot_dat = plot_dat %>% dplyr::rename("Beta" = "Estimate", "StEr" = "Std..Error", "p" = "Pr...t..") %>% dplyr::select(Beta, StEr, Variables, Metric,p)

beta_plot = ggplot(plot_dat, aes(x = Variables, y = Beta, group = Metric)) + 
  geom_point(position = position_dodge(width = 0.4), aes(color= Metric), size = 1.8)+
  geom_errorbar(aes(ymin=Beta-StEr, ymax=Beta+StEr, y=Beta, x = (Variables), color = Metric), position = position_dodge(width = 0.4), width = 0.1)+
  labs(x="", y = "Beta ± Standard Error")+
  #guides(color = guide_legend(reverse = FALSE))+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 45))+
  ggtitle("Predictors of Brain Age controlling for Age, Sex, and the Age-Sex Interaction")

# we plot the interaction effect for the multimodal model (significant interaction)
combi_int = allEffects(combi_model)
plot(combi_int)


# 
# 
# df2.2 = df2
# levels(df2.2$Handedness) = c("right","left","ambi")
# # make violin plots for age and predicted age for T1 data
# 
# # Basic violin plot T1 predicted age
# age_plot = ggplot(na.omit(df2.2), aes(age, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="Chronological Age",x="Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal() +theme(legend.position = "none")
# # Basic violin plot T1 predicted age
# L_T1 = ggplot(na.omit(df2.2), aes(LT1, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="L Hemisphere T1w",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# R_T1 = ggplot(na.omit(df2.2), aes(RT1, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="R Hemisphere T1w",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# T1_plot = ggarrange(L_T1,R_T1, common.legend = T, legend = "none")
# # Basic violin plot dMRI predicted age
# L_dMRI = ggplot(na.omit(df2.2), aes(LdMRI, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="L Hemisphere dMRI",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# R_dMRI = ggplot(na.omit(df2.2), aes(RdMRI, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="R Hemisphere dMRI",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# dMRI_plot = ggarrange(L_dMRI,R_dMRI, common.legend = T, legend = "none")
# # Basic violin plot combi predicted age
# L_combi = ggplot(na.omit(df2.2), aes(Lcombi, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="L Hemi multimodal",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# R_combi = ggplot(na.omit(df2.2), aes(Rcombi, Handedness, fill = Handedness)) + 
#   geom_violin(trim=FALSE)+
#   geom_boxplot(width=0.1, position = position_dodge(0.5))+
#   labs(title="R Hemi multimodal",x="Predicted Age", y = "Handedness") + 
#   scale_fill_brewer(palette="Dark2") + theme_minimal()
# multimodal_plot = ggarrange(L_combi,R_combi, common.legend = T, legend = "none")
# 
# violins = ggarrange(age_plot, T1_plot,dMRI_plot,multimodal_plot)
# Plot2 = ggarrange(beta_plot,violins)
# 
# ggsave("/tsd/p33/home/p33-maxk/export/raincloud.pdf",Plot2, width = 18, height = 10)
# ggsave("/cluster/projects/p33/users/maxk/UKB/export/Hemi_Plot2.pdf",Plot2, width = 18, height = 10)


# useful for common axes:
# require(grid)   # for the textGrob() function
# 
# figure <- ggarrange(L_T1 + rremove("ylab") + rremove("xlab"), R_T1 + rremove("ylab") + rremove("xlab"), #gg_3 + rremove("ylab") + rremove("xlab"), gg_4+ rremove("ylab") + rremove("xlab"), # remove axis labels from plots
#                     labels = NULL,
#                     ncol = 2, nrow = 2,
#                     common.legend = TRUE, legend = "bottom",
#                     align = "hv", 
#                     font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
# 
# annotate_figure(figure, left = textGrob("Common y-axis", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
#                 bottom = textGrob("Common x-axis", gp = gpar(cex = 1.3)))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 

## HEMISPHERIC DIFFERENCES ####

########### PREP
# check whether there are differences when holding sex, age, and scanner site constant
covars = demo %>% dplyr::select(eid, sex, site)
df.1 = merge(covars, df, by = "eid")

# BAG differences
## melt data frames into long format
df3 = df.1 %>% dplyr::select(-eid, -age, -sex, -site)
df2 = df.1 %>% dplyr::select(-eid, -sex, -site)
df2 = melt(df2)
df3 = melt(df3)

## add covariates to melted data frames
df3$sex = c(replicate(9, df.1$sex))
df3$site = c(replicate(9, df.1$site))
df3$age = c(replicate(9, df.1$age))

######## DESCRIBE

## describe
df2 %>% group_by(variable) %>% summarise(across(.fns = c(mean, sd)))
df3 %>% group_by(variable) %>% summarise(across(.fns = c(mean, sd)))

# these differences are significant (across all brain age estimates)
kruskal.test(value~variable, data = df2)
kruskal.test(value~variable, data = df3)

# Crude BAG differences are not significant when calculating hemispheric BAG from single modalities
t.test(df$LdMRI,df$RdMRI)
t.test(df$LT1,df$RT1)
# predicted age differences
t.test(df0$LT1,df0$RT1)
t.test(df0$LdMRI,df0$RdMRI)
t.test(df0$Lcombi,df0$Rcombi)

# but differences across modalities within each hemisphere are significant
t.test(df$LdMRI, df$LT1)
t.test(df$RdMRI, df$RT1)
# and the hemispheric differences in BAG differ between T1 and dwMRI
t.test(df$LdMRI-df$RdMRI, df$LT1-df$RT1)

Ldf = df0 %>% dplyr::select(LT1, LdMRI, Lcombi)
Rdf = df0 %>% dplyr::select(RT1, RdMRI, Rcombi)
Ldf = melt(Ldf)
Rdf = melt(Rdf)
Ldf_mod = aov(value ~ variable, Ldf)
summary(Ldf_mod)
kruskal.test(value ~ variable, Ldf) # check whether Kruskal Wallis yields teh same result

Rdf_mod = aov(value ~ variable, Rdf)
summary(Rdf_mod)
kruskal.test(value ~ variable, Rdf) # check whether Kruskal Wallis yields teh same result

# post hoc
Rdf %>% pairwise_wilcox_test(value ~ variable)
Rdf %>% wilcox_effsize(value ~ variable)
Ldf %>% pairwise_wilcox_test(value ~ variable)
Ldf %>% wilcox_effsize(value ~ variable)

# plot (data prep)
df2 = df
df2$T1 = (df2$LT1-df2$RT1)
df2$dwMRI = (df2$LdMRI-df2$RdMRI)
df3 = df2[,6:7]
df3 = melt(df3)

# density plots for hemispheric BAG differences
ggplot(data=df3, aes(x=value, group=variable, fill=variable)) +
  geom_density(adjust=1.5, alpha=.4) #+ theme_ipsum()

# hemispheric differences in BAG are not strongly correlated (to each other or single hemisphere BAGs)
df2 %>% dplyr::select(-eid)%>%cor()

df4 = merge(df2, demo, by = "eid")

# run linear mixed models with age, sex (0=f, 1=m) and scanner site
model1 = lmer(LdMRI ~ sex*age + (1|site), data = df4)
model2 = lmer(RdMRI ~ sex*age + (1|site), data = df4)
summary(model1)
summary(model2)
model3 = lmer(LT1 ~ sex*age + (1|site), data = df4)
model4 = lmer(RT1 ~ sex*age + (1|site), data = df4)
summary(model3)
summary(model4)

# age, sex and assessment centre explain more of the single hemisphere BAG in dwMRI
r.squaredGLMM(model1)
r.squaredGLMM(model2)
r.squaredGLMM(model3)
r.squaredGLMM(model4)

# sex effects are the opposite in T1 vs dwMRI: females white matter is later affected than in males, for grey matter it is the opposite
t.test(df4$LdMRI ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(LdMRI))
t.test(df4$RdMRI ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(RdMRI))
t.test(df4$LT1 ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(RT1))
t.test(df4$RT1 ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(LT1))

# also hemispheric BAG differences effect sizes are bigger for dwMRI
cohens_d(df4$LdMRI ~ df4$sex, alternative = "two.sided")
cohens_d(df4$RdMRI ~ df4$sex, alternative = "two.sided")
cohens_d(df4$LT1 ~ df4$sex, alternative = "two.sided")
cohens_d(df4$RT1 ~ df4$sex, alternative = "two.sided")

# what about their difference? >> then it's the same direction ...
t.test(df4$LdMRI-df4$RdMRI ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(LdMRI-RdMRI))
t.test(df4$LT1-df4$RT1 ~ df4$sex, alternative = "two.sided")
df4 %>% group_by(sex) %>% summarize(SD = sd(LT1-RT1))
# and effect sizes are more or less the same...
cohens_d(df4$LdMRI-df4$RdMRI ~ df4$sex, alternative = "two.sided")
cohens_d(df4$LT1-df4$RT1 ~ df4$sex, alternative = "two.sided")


# relationships of hemispheric brain age with health outcomes
df5 = dplyr::rename(df4, BAG = LT1)
df6 = dplyr::rename(df4, BAG = RT1)
df7 = dplyr::rename(df4, BAG = LdMRI)
df8 = dplyr::rename(df4, BAG = RdMRI)
df9= dplyr::rename(df4, BAG = T1)
df10= dplyr::rename(df4, BAG = dwMRI)
pred_list = list(df5,df6,df7,df8,df9,df10)
health = list()
cognition = list()
#for (i in 1:length(pred_list)){
#  health[[i]] = lmer(BAG~sex*age+WHR+birth_weight+diabetic+hypertension+(1|site), pred_list[[i]])
#  cognition[[i]] = lmer(BAG~sex*age+prospective_memory+fluid_intelligence+matrix_puzzles_solved+digits_remembered+digit_sub_correct+(1|site), pred_list[[i]])
#}
