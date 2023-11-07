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
library(emmeans)
#install.packages("dataPreparation")
library(dataPreparation)

# read in data
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/T1.csv")
dMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/dMRI.csv")
#multimodal = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/combi.csv")
demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")
# HEMISHPERIC DIFFERENCES ####
# select the data for each hemisphere and get mean scores for each of the metrics
# T1: surface area, thickness, volume
T1_lh = T1 %>% dplyr::select(starts_with(c("lh_")))
T1_rh = T1 %>% dplyr::select(starts_with(c("rh_")))

T1_t.tests = list()
for (i in 1:ncol(T1_lh)){
  x = T1_lh[i]
  y = T1_rh[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  T1_t.tests[[i]] = t.test(df$x, df$y, paired = T)
}

T1_p = c()
for (i in 1:length(T1_t.tests)){
  T1_p[i] = (T1_t.tests[[i]]$p.value)
}
# absolute value
sum(ifelse(T1_p > 0.05/length(T1_t.tests), 1, 0))
# relative value
sum(ifelse(T1_p > 0.05/length(T1_t.tests), 1, 0))/length(T1_t.tests)


# dMRI: 28 metrics
dMRI_lh = dMRI %>% dplyr::select(ends_with(c("L")))
dMRI_rh = dMRI %>% dplyr::select(ends_with(c("R")))

dMRI_t.tests = list()
for (i in 1:ncol(dMRI_lh)){
  x = dMRI_lh[i]
  y = dMRI_rh[i]
  df = data.frame(x,y)
  colnames(df) = c("x","y")
  dMRI_t.tests[[i]] = t.test(df$x, df$y, paired = T)
}

dMRI_p = c()
for (i in 1:length(dMRI_t.tests)){
  dMRI_p[i] = (dMRI_t.tests[[i]]$p.value)
}
# absolute value
sum(ifelse(dMRI_p > 0.05/length(dMRI_t.tests), 1, 0))
# relative value
sum(ifelse(dMRI_p > 0.05/length(dMRI_t.tests), 1, 0))/length(dMRI_t.tests)


# AGE DEPENDENCIES OF MEAN METRICS ####

#T1w
T1_lh_thickness = T1_lh %>% dplyr::select(contains(c("thickness"))) %>% rowMeans()
T1_lh_volume = T1_lh %>% dplyr::select(contains(c("volume"))) %>% rowMeans()
T1_lh_area = T1_lh %>% dplyr::select(contains(c("area"))) %>% rowMeans()
T1_rh_thickness = T1_rh %>% dplyr::select(contains(c("thickness"))) %>% rowMeans()
T1_rh_volume = T1_rh %>% dplyr::select(contains(c("volume"))) %>% rowMeans()
T1_rh_area = T1_rh %>% dplyr::select(contains(c("area"))) %>% rowMeans()
T1_df = data.frame(T1$eid, T1_lh_thickness,T1_rh_thickness, T1_lh_area, T1_rh_area, T1_lh_volume, T1_rh_volume)

# dMRI L
BRIA_vintra_lh = dMRI_lh %>% dplyr::select(contains(c("v_intra_"))) %>% rowMeans()
BRIA_vextra_lh = dMRI_lh %>% dplyr::select(contains(c("v_extra_"))) %>% rowMeans()
BRIA_vcsf_lh = dMRI_lh %>% dplyr::select(contains(c("v_csf_"))) %>% rowMeans()
BRIA_micrord_lh = dMRI_lh %>% dplyr::select(contains(c("micro_Rd_"))) %>% rowMeans()
BRIA_microfa_lh = dMRI_lh %>% dplyr::select(contains(c("micro_FA_"))) %>% rowMeans()
BRIA_microax_lh = dMRI_lh %>% dplyr::select(contains(c("micro_Ax_"))) %>% rowMeans()
BRIA_microadc_lh = dMRI_lh %>% dplyr::select(contains(c("micro_ADC_"))) %>% rowMeans()
BRIA_dradextra_lh = dMRI_lh %>% dplyr::select(contains(c("Drad_extra"))) %>% rowMeans()
BRIA_daxintra_lh = dMRI_lh %>% dplyr::select(contains(c("Dax_intra"))) %>% rowMeans()
BRIA_daxextra_lh = dMRI_lh %>% dplyr::select(contains(c("Dax_extra_"))) %>% rowMeans()
DTI_AD_lh = dMRI_lh %>% dplyr::select(contains(c("AD_"))) %>% rowMeans()
DTI_RD_lh = dMRI_lh %>% dplyr::select(contains(c("RD_"))) %>% rowMeans()
DTI_MD_lh = dMRI_lh %>% dplyr::select(contains(c("MD_"))) %>% rowMeans()
DTI_FA_lh = dMRI_lh %>% dplyr::select(contains(c("FA_"))) %>% rowMeans()

DKI_AK_lh = dMRI_lh %>% dplyr::select(contains(c("ak_"))) %>% rowMeans()
DKI_RK_lh = dMRI_lh %>% dplyr::select(contains(c("rk_"))) %>% rowMeans()
DKI_MK_lh = dMRI_lh %>% dplyr::select(contains(c("mk_"))) %>% rowMeans()

SMT_FA_lh = dMRI_lh %>% dplyr::select(contains(c("smt_fa_"))) %>% rowMeans()
SMT_MD_lh = dMRI_lh %>% dplyr::select(contains(c("smt_md_"))) %>% rowMeans()
SMT_trans_lh = dMRI_lh %>% dplyr::select(contains(c("smt_trans_"))) %>% rowMeans()
SMT_long_lh = dMRI_lh %>% dplyr::select(contains(c("smt_long"))) %>% rowMeans()

SMTmc_d_lh  = dMRI_lh %>% dplyr::select(contains(c("smt_mc_diff"))) %>% rowMeans()
SMTmc_extramd_lh = dMRI_lh %>% dplyr::select(contains(c("smt_mc_extramd"))) %>% rowMeans()
SMTmc_extratrans_lh = dMRI_lh %>% dplyr::select(contains(c("smt_mc_extratrans"))) %>% rowMeans()
SMTmc_intra_lh = dMRI_lh %>% dplyr::select(contains(c("smt_mc_intra"))) %>% rowMeans()

WMTI_awf_lh= dMRI_lh %>% dplyr::select(contains(c("awf"))) %>% rowMeans()
WMTI_radead_lh= dMRI_lh %>% dplyr::select(contains(c("radEAD"))) %>% rowMeans()
WMTI_axead_lh= dMRI_lh %>% dplyr::select(contains(c("axEAD"))) %>% rowMeans()

# dMRI R
BRIA_vintra_rh = dMRI_rh %>% dplyr::select(contains(c("v_intra_"))) %>% rowMeans()
BRIA_vextra_rh = dMRI_rh %>% dplyr::select(contains(c("v_extra_"))) %>% rowMeans()
BRIA_vcsf_rh = dMRI_rh %>% dplyr::select(contains(c("v_csf_"))) %>% rowMeans()
BRIA_micrord_rh = dMRI_rh %>% dplyr::select(contains(c("micro_Rd_"))) %>% rowMeans()
BRIA_microfa_rh = dMRI_rh %>% dplyr::select(contains(c("micro_FA_"))) %>% rowMeans()
BRIA_microax_rh = dMRI_rh %>% dplyr::select(contains(c("micro_Ax_"))) %>% rowMeans()
BRIA_microadc_rh = dMRI_rh %>% dplyr::select(contains(c("micro_ADC_"))) %>% rowMeans()
BRIA_dradextra_rh = dMRI_rh %>% dplyr::select(contains(c("Drad_extra"))) %>% rowMeans()
BRIA_daxintra_rh = dMRI_rh %>% dplyr::select(contains(c("Dax_intra"))) %>% rowMeans()
BRIA_daxextra_rh = dMRI_rh %>% dplyr::select(contains(c("Dax_extra_"))) %>% rowMeans()

DTI_AD_rh = dMRI_rh %>% dplyr::select(contains(c("AD_"))) %>% rowMeans()
DTI_RD_rh = dMRI_rh %>% dplyr::select(contains(c("RD_"))) %>% rowMeans()
DTI_MD_rh = dMRI_rh %>% dplyr::select(contains(c("MD_"))) %>% rowMeans()
DTI_FA_rh = dMRI_rh %>% dplyr::select(contains(c("FA_"))) %>% rowMeans()

DKI_AK_rh = dMRI_rh %>% dplyr::select(contains(c("ak_"))) %>% rowMeans()
DKI_RK_rh = dMRI_rh %>% dplyr::select(contains(c("rk_"))) %>% rowMeans()
DKI_MK_rh = dMRI_rh %>% dplyr::select(contains(c("mk_"))) %>% rowMeans()

SMT_FA_rh = dMRI_rh %>% dplyr::select(contains(c("smt_fa_"))) %>% rowMeans()
SMT_MD_rh = dMRI_rh %>% dplyr::select(contains(c("smt_md_"))) %>% rowMeans()
SMT_trans_rh = dMRI_rh %>% dplyr::select(contains(c("smt_trans_"))) %>% rowMeans()
SMT_long_rh = dMRI_rh %>% dplyr::select(contains(c("smt_long"))) %>% rowMeans()

SMTmc_d_rh  = dMRI_rh %>% dplyr::select(contains(c("smt_mc_diff"))) %>% rowMeans()
SMTmc_extramd_rh = dMRI_rh %>% dplyr::select(contains(c("smt_mc_extramd"))) %>% rowMeans()
SMTmc_extratrans_rh = dMRI_rh %>% dplyr::select(contains(c("smt_mc_extratrans"))) %>% rowMeans()
SMTmc_intra_rh = dMRI_rh %>% dplyr::select(contains(c("smt_mc_intra"))) %>% rowMeans()

WMTI_awf_rh= dMRI_rh %>% dplyr::select(contains(c("awf"))) %>% rowMeans()
WMTI_radead_rh= dMRI_rh %>% dplyr::select(contains(c("radEAD"))) %>% rowMeans()
WMTI_axead_rh= dMRI_rh %>% dplyr::select(contains(c("axEAD"))) %>% rowMeans()

# merge all data
dMRI_df = data.frame(BRIA_vintra_lh,BRIA_vintra_rh,BRIA_vextra_lh,BRIA_vextra_rh,
                     BRIA_vcsf_lh,BRIA_vcsf_rh,BRIA_micrord_lh,BRIA_micrord_rh,
                     BRIA_microfa_lh,BRIA_microfa_rh,
                     BRIA_microax_lh,BRIA_microax_rh,
                     BRIA_microadc_lh,BRIA_microadc_rh,
                     BRIA_dradextra_lh,BRIA_dradextra_rh,
                     BRIA_daxintra_lh, BRIA_daxintra_rh,
                     BRIA_daxextra_lh,  BRIA_daxextra_rh, 
                     DKI_AK_lh, DKI_AK_rh, 
                     DKI_RK_lh, DKI_RK_rh, 
                     DKI_MK_lh,DKI_MK_rh,
                     DTI_AD_lh, DTI_AD_rh, 
                     DTI_RD_lh,  DTI_RD_rh, 
                     DTI_MD_lh, DTI_MD_rh, 
                     DTI_FA_lh,DTI_FA_rh,
                     SMT_FA_lh,SMT_FA_rh,
                     SMT_MD_lh,SMT_MD_rh,
                     SMT_trans_lh,SMT_trans_rh,
                     SMT_long_lh,SMT_long_rh,
                     SMTmc_d_lh, SMTmc_d_rh,
                     SMTmc_extramd_lh,SMTmc_extramd_rh,
                     SMTmc_extratrans_lh,SMTmc_extratrans_rh,
                     SMTmc_intra_lh,SMTmc_intra_rh,
                     WMTI_awf_lh,WMTI_awf_rh,
                     WMTI_radead_lh,WMTI_radead_rh,
                     WMTI_axead_lh,WMTI_axead_rh, dMRI$eid)
names(dMRI_df) = c("BRIA_vintra_lh","BRIA_vintra_rh","BRIA_vextra_lh","BRIA_vextra_rh",
                   "BRIA_vcsf_lh","BRIA_vcsf_rh","BRIA_micrord_lh","BRIA_micrord_rh",
                   "BRIA_microfa_lh","BRIA_microfa_rh",
                   "BRIA_microax_lh","BRIA_microax_rh",
                   "BRIA_microadc_lh","BRIA_microadc_rh",
                   "BRIA_dradextra_lh","BRIA_dradextra_rh",
                   "BRIA_daxintra_lh", "BRIA_daxintra_rh",
                   "BRIA_daxextra_lh",  "BRIA_daxextra_rh", 
                   "DKI_AK_lh", "DKI_AK_rh", 
                   "DKI_RK_lh", "DKI_RK_rh", 
                   "DKI_MK_lh","DKI_MK_rh",
                   "DTI_AD_lh", "DTI_AD_rh", 
                   "DTI_RD_lh",  "DTI_RD_rh", 
                   "DTI_MD_lh", "DTI_MD_rh", 
                   "DTI_FA_lh","DTI_FA_rh",
                   "SMT_FA_lh","SMT_FA_rh",
                   "SMT_MD_lh","SMT_MD_rh",
                   "SMT_trans_lh","SMT_trans_rh",
                   "SMT_long_lh","SMT_long_rh",
                   "SMTmc_d_lh", "SMTmc_d_rh",
                   "SMTmc_extramd_lh","SMTmc_extramd_rh",
                   "SMTmc_extratrans_lh","SMTmc_extratrans_rh",
                   "SMTmc_intra_lh","SMTmc_intra_rh",
                   "WMTI_awf_lh","WMTI_awf_rh",
                   "WMTI_radead_lh","WMTI_radead_rh",
                   "WMTI_axead_lh","WMTI_axead_rh", "eid"
)

dMRI_df = merge(demo,dMRI_df, by = "eid")
tab = dMRI_df %>% dplyr::select(-eid, -Sex, -Scanner)%>%na.omit() %>%cor()

# make a single df for the plots
T1_df = T1_df %>% dplyr::rename("eid" = "T1.eid") 
orig_df = merge(dMRI_df,  T1_df, by = "eid")

# remove outliers that are 3 SDs over mean
# this is an important step as otherwise the whole curve-fitting results in nonsense

plot_df = remove_sd_outlier(orig_df)

# no, we can correct these metrics for sex, age, and scanner site
# for that, first list all the metrics
metrics_list = c("BRIA_vintra_lh","BRIA_vintra_rh","BRIA_vextra_lh","BRIA_vextra_rh",
                 "BRIA_vcsf_lh","BRIA_vcsf_rh","BRIA_micrord_lh","BRIA_micrord_rh",
                 "BRIA_microfa_lh","BRIA_microfa_rh",
                 "BRIA_microax_lh","BRIA_microax_rh",
                 "BRIA_microadc_lh","BRIA_microadc_rh",
                 "BRIA_dradextra_lh","BRIA_dradextra_rh",
                 "BRIA_daxintra_lh", "BRIA_daxintra_rh",
                 "BRIA_daxextra_lh",  "BRIA_daxextra_rh", 
                 "DKI_AK_lh", "DKI_AK_rh", 
                 "DKI_RK_lh", "DKI_RK_rh", 
                 "DKI_MK_lh","DKI_MK_rh",
                 "DTI_AD_lh", "DTI_AD_rh", 
                 "DTI_RD_lh",  "DTI_RD_rh", 
                 "DTI_MD_lh", "DTI_MD_rh", 
                 "DTI_FA_lh","DTI_FA_rh",
                 "SMT_FA_lh","SMT_FA_rh",
                 "SMT_MD_lh","SMT_MD_rh",
                 "SMT_trans_lh","SMT_trans_rh",
                 "SMT_long_lh","SMT_long_rh",
                 "SMTmc_d_lh", "SMTmc_d_rh",
                 "SMTmc_extramd_lh","SMTmc_extramd_rh",
                 "SMTmc_extratrans_lh","SMTmc_extratrans_rh",
                 "SMTmc_intra_lh","SMTmc_intra_rh",
                 "WMTI_awf_lh","WMTI_awf_rh",
                 "WMTI_radead_lh","WMTI_radead_rh",
                 "WMTI_axead_lh","WMTI_axead_rh",
                 "T1_lh_thickness","T1_rh_thickness",
                 "T1_lh_area", "T1_rh_area",
                 "T1_lh_volume", "T1_rh_volume")

prediction_df = data.frame(matrix(ncol = length(metrics_list), nrow = nrow(plot_df)))
mods = list()
#betaMean = c()
#psMean = c()
for (i in 1:length(metrics_list)){
  formula1 = paste(metrics_list[i], "~ Age+Sex+(1|Scanner)")
  mods[[i]] = lmer(formula1, data = plot_df)
  #betaMean[i] = summary(mods[[i]])$coefficients[2,1]
  #psMean[i] = summary(mods[[i]])$coefficients[2,5]*length(metrics_list)
  prediction_df[i] = predict(mods[[i]])
}

# get abs LI
LEFT = plot_df %>% dplyr::select(contains("lh"))
RIGHT = plot_df %>% dplyr::select(contains("rh"))
LI = function(L,R){
  abs((L-R)/(L+R))
}
LIdata = data.frame(matrix(ncol = ncol(LEFT),nrow = nrow(LEFT)))
for (i in 1:ncol(LEFT)){
  LIdata[i] = scale(LI(LEFT[[i]], RIGHT[[i]]))
}
colnames(LIdata) = colnames(LEFT)
LIdata = LIdata %>% rename_with(~str_remove(., '_lh'))
LIdata$Age = plot_df$Age
LIdata$Sex = plot_df$Sex
LIdata$Scanner = plot_df$Scanner

# now loop over it
mods = list()
betaMean = c()
psMean = c()
for (i in 1:ncol(LEFT)){
  formula1 = paste(colnames(LIdata)[i], "~ Age+Sex+(1|Scanner)")
  mods[[i]] = lmer(formula1, data = LIdata)
  betaMean[i] = summary(mods[[i]])$coefficients[2,1]
  psMean[i] = summary(mods[[i]])$coefficients[2,5]*ncol(LEFT)
  #prediction_df[i] = predict(mods[[i]])
}
# get model metrics
padjMean = p.adjust(psMean, method = "bonferroni")
mean_LI = data.frame(colnames(LIdata[1:31]),betaMean,padjMean)
mean(abs(mean_LI$betaMean))
sd(abs(mean_LI$betaMean))

# add column names for the metrics and age for plotting
colnames(prediction_df)=metrics_list
# now, standardize values
prediction_df=prediction_df%>%mutate_if(is.numeric,scale)

# add age to prediction df
prediction_df$age = plot_df$Age


# now we can plot the corrected values
plot_list = list()
for (i in 1:length(metrics_list)){
  plot1_df = data.frame(prediction_df$age, prediction_df[i])
  colnames(plot1_df) = c("age", "pred")
  plot_list[[i]] =  ggplot(plot1_df, aes(x = age, y = pred)) +
    #geom_bin2d(bins = 150)+
    #scale_fill_continuous(type = "viridis")+
    stat_smooth(method = "gam", col = "red") +
    labs(x = "Age", y = paste(metrics_list[i]))+
    theme_bw()+ theme(legend.position = "none")
}
# arrange the plots as 2x31 & save the resulting plot
age_curves = ggarrange(plotlist = plot_list) %>%
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_Supplement7.pdf", plot = ., width = 15, height = 25, device = "pdf")
  #ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_age_curves_LMER.pdf", plot = ., width = 15, height = 25, device = "pdf")

# save also raw plots
prediction_df2 = plot_df %>% dplyr::select(-eid,-Age,-Sex,-Scanner)

# now we plot the UNcorrected values WITHOUT SCATTER
plot_list2 = list()
for (i in 1:(length(metrics_list))){
  temp_df = data.frame(pred = prediction_df2[[i]],age = plot_df$Age)
  yname=metrics_list[i]
  #names(temp_df) = c("age", "pred")
  plot_list2[[i]] =  ggplot(temp_df, aes(age,pred)) +
  #  geom_bin2d(bins = 150)+
  #  scale_fill_continuous(type = "viridis")+
    stat_smooth(method = "gam", col = "blue") +
    labs(x = "Age", y = yname)+
    theme_bw()+ theme(legend.position = "none")
}

# arrange the plots as 2x31 & save the resulting plot
age_curves = ggarrange(plotlist = plot_list2) %>%
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_Supplement6.pdf", plot = ., width = 15, height = 25, device = "pdf")
#ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_age_curves.pdf", plot = ., width = 15, height = 25, device = "pdf")

# WITH SCATTER
plot_list2 = list()
for (i in 1:(length(metrics_list))){
  temp_df = data.frame(pred = prediction_df2[[i]],age = plot_df$Age)
  yname=metrics_list[i]
  #names(temp_df) = c("age", "pred")
  plot_list2[[i]] =  ggplot(temp_df, aes(age,pred)) +
    geom_bin2d(bins = 150)+
    scale_fill_continuous(type = "viridis")+
    stat_smooth(method = "gam", col = "blue") +
    labs(x = "Age", y = yname)+
    theme_bw()+ theme(legend.position = "none")
}

# plot the uncorrected relationships
age_curves = ggarrange(plotlist = plot_list2) %>% 
  #ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_age_curves.pdf", plot = ., width = 15, height = 25, device = "pdf")
  ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_Supplement14.pdf", plot = ., width = 15, height = 25, device = "pdf")

## IN-TEXT PLOT
# prep
LEFT = prediction_df %>% dplyr::select(contains("lh"), age)
RIGHT = prediction_df %>% dplyr::select(contains("rh"), age)
names(LEFT) = c("BRIA_vintra",      "BRIA_vextra"  ,    "BRIA_vcsf"     ,   "BRIA_micrord"   ,  "BRIA_microfa"   ,  "BRIA_microax"  ,   "BRIA_microadc"  , 
                "BRIA_dradextra" ,  "BRIA_daxintra"  ,  "BRIA_daxextra"  ,  "DKI_AK"     ,      "DKI_RK"     ,      "DKI_MK"     ,      "DTI_AD"  ,        
                "DTI_RD"    ,       "DTI_MD"    ,       "DTI_FA"   ,        "SMT_FA"     ,      "SMT_MD"       ,    "SMT_trans"    ,    "SMT_long"   ,     
                "SMTmc_d"      ,    "SMTmc_extramd" ,   "SMTmc_extratrans", "SMTmc_intra"  ,    "WMTI_awf"     ,    "WMTI_radead"    ,  "WMTI_axead"      ,
                "T1_thickness"   ,  "T1_area"  ,        "T1_volume"     ,   "age")
names(RIGHT) = names(LEFT)
both = rbind(LEFT, RIGHT)
both$Hemisphere = c(replicate(nrow(LEFT),"L"),replicate(nrow(RIGHT),"R"))
betternames = c("BRIA - V intra",      "BRIA - V extra"  ,    "BRIA - V CSF"     ,   "BRIA - microRD"   ,  "BRIA - microFA"   ,  "BRIA - microAX"  ,   "BRIA - microADC"  , 
                              "BRIA - dradextra" ,  "BRIA - daxintra"  ,  "BRIA - daxextra"  ,  "DKI - AK"     ,      "DKI - RK"     ,      "DKI - MK"     ,      "DTI - AD"  ,        
                              "DTI - RD"    ,       "DTI - MD"    ,       "DTI - FA"   ,        "SMT - FA"     ,      "SMT - MD"       ,    "SMT - trans"    ,    "SMT - long"   ,     
                              "SMTmc - diff"      ,    "SMTmc - extraMD" ,   "SMTmc - extra trans", "SMTmc - intra"  ,    "WMTI - AWF"     ,    "WMTI - radEAD"    ,  "WMTI - axEAD"      ,
                              "T1 - thickness"   ,  "T1 - area"  ,        "T1 - volume")
# plot
plotlist3=list()
for (i in 1:31){
  temp_df = data.frame(pred = both[[i]],age = both$age, Hemisphere = both$Hemisphere)
  yname=betternames[i]
  #names(temp_df) = c("age", "pred")
  plotlist3[[i]] =  ggplot() +
  #  geom_bin2d(bins = 150)+
  #  scale_fill_continuous(type = "viridis")+
  stat_smooth(data = temp_df, aes(x=age,y=pred, group = Hemisphere, color = Hemisphere),method = "gam") +
  labs(x = "Age", y = yname)+
  theme_bw()+ theme(legend.position = "none")
}


age_curves = ggarrange(plotlist = plotlist3, common.legend = T) %>% 
  #ggsave(filename = "/tsd/p33/home/p33-maxk/export/Hemi_age_curves.pdf", plot = ., width = 15, height = 25, device = "pdf")
ggsave(filename = "/cluster/projects/p33/users/maxk/UKB/export/Hemi_Fig1.pdf", plot = ., width = 15, height = 25, device = "pdf")



# MODEL AGE SENSITIVITY ####
# run models predicting age with and without mean hemispheric grey and white matter metrics
mod0list = list()
mod1list = list()
for (i in (metrics_list)){
  f0 = formula(paste("age~sex+site"))
  f1 = formula(paste("age~",i, "+sex+site"))
  mod0list[[]] = lm(f0,data = plot_df)
  mod1list = lm(f1,data = plot_df)
}
# SAMPLE DESCRIPTOR ####

# multimodal
multimodal = merge(multimodal, dMRI, by = "eid")
prop.table(table(multimodal$site))
prop.table(table(multimodal$sex))
max(na.omit(multimodal$age.x))
min(multimodal$age.x)
mean(multimodal$age.x)
sd(multimodal$age.x)
median(multimodal$age.x)
mad(multimodal$age.x)

# T1w
prop.table(table(orig_df$site))
prop.table(table(orig_df$sex))
max(orig_df$age)
min(orig_df$age)
mean(orig_df$age)
sd(orig_df$age)
median(orig_df$age)
mad(orig_df$age)

# dMRI
prop.table(table(dMRI$site))
prop.table(table(dMRI$sex))
max(dMRI$age)
min(dMRI$age)
mean(dMRI$age)
sd(dMRI$age)
median(dMRI$age)
mad(dMRI$age)
