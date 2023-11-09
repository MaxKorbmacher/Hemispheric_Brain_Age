# Hemispheric brain age plot: beta coefficients by p-values EXTENSION: ABSOLUTE ASYMMETRIES
# Max Korbmacher, 05.09.2023, R version 4.1.2 (2021-11-01), updated 09.11.2023
####### PREP ####

# load library
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, stringr, ggpubr, ggrepel, Rmpfr, egg)

#set working directory
setwd("/home/max/Downloads/")
# load data (both sexes and sex-stratified)
T1 = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_both_sexes.csv")
dmri =  read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_both_sexes.csv")
T1_m = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_MALES.csv")
T1_f= read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_T1_FEMALES.csv")
dmri_m = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_MALES.csv")
dmri_f = read.csv("Hemi_NEW_REGIONAL_LINEAR_LI_age_ass_dMRI_FEMALES.csv")
colnames(T1) = c("names", "feature", "beta_out", "t")
colnames(T1_m) = colnames(T1)
colnames(T1_f) = colnames(T1)
colnames(dmri) = colnames(T1)
colnames(dmri_m) = colnames(T1)
colnames(dmri_f) = colnames(T1)

# estimate and adjust p-values
## make function for precise p-val display (no rounding) as estimated from t-values
.N <- function(.) mpfr(., precBits = 200)
# T1w
T1$pvals = 2 * pnorm(-abs(.N(T1$t)))
T1$p_adj = formatMpfr(mpfr((T1$pvals)*nrow(T1),200))
T1_m$pvals = 2 * pnorm(-abs(.N(T1_m$t)))
T1_m$p_adj = mpfr(T1_m$pvals*nrow(T1_m),200)
T1_f$pvals = 2 * pnorm(-abs(.N(T1_f$t)))
T1_f$p_adj = mpfr(T1_f$pvals*nrow(T1_f),200)
# dMRI
dmri$pvals = 2 * pnorm(-abs(.N(dmri$t)))
dmri$p_adj = mpfr(dmri$pvals*nrow(dmri),200)
dmri_m$pvals = 2 * pnorm(-abs(.N(dmri_m$t)))
dmri_m$p_adj = mpfr(dmri_m$pvals*nrow(dmri_m),200)
dmri_f$pvals = 2 * pnorm(-abs(.N(dmri_f$t)))
dmri_f$p_adj = mpfr(dmri_f$pvals*nrow(dmri_f),200)
# limit adjusted values to range 0 to 1 (not bigger than 1)
T1$p_adj=ifelse(T1$p_adj>1,1,T1$p_adj)
T1_m$p_adj=ifelse(T1_m$p_adj>1,1,T1_m$p_adj)
T1_f$p_adj=ifelse(T1_f$p_adj>1,1,T1_f$p_adj)
dmri$p_adj=ifelse(dmri$p_adj>1,1,dmri$p_adj)
dmri_m$p_adj=ifelse(dmri_m$p_adj>1,1,dmri_m$p_adj)
dmri_f$p_adj=ifelse(dmri_f$p_adj>1,1,dmri_f$p_adj)

###### STATS ####
# dmri %>% filter((p_adj) < 0.05) %>% summarize(Mean = mean(abs(beta_out)), SD = sd(abs(beta_out)))
# T1 %>% filter(p_adj < 0.05) %>% summarize(Mean = mean(abs(beta_out)), SD = sd(abs(beta_out)))
# rbind(dmri,T1) %>% filter(p_adj < 0.05) %>% summarize(Mean = mean(abs(beta_out)), SD = sd(abs(beta_out)))

####### PLOT PREP (both sexes) ####
# log10 transform adjusted p_values
T1$log10p = log10(2 * pnorm(-abs(.N(T1$t))))*-1
T1_f$log10p = log10(2 * pnorm(-abs(.N(T1_f$t))))*-1
T1_m$log10p = log10(2 * pnorm(-abs(.N(T1_m$t))))*-1
dmri$log10p = log10(2 * pnorm(-abs(.N(dmri$t))))*-1
dmri_m$log10p = log10(2 * pnorm(-abs(.N(dmri_m$t))))*-1
dmri_f$log10p = log10(2 * pnorm(-abs(.N(dmri_f$t))))*-1

# remove prefixes and suffixes (remnant)
## dMRI features
dmri$feature=str_sub(dmri$feature,1,-2)
dmri$feature=gsub("BRIA.BRIA.","BRIA.",dmri$feature)
dmri$feature=gsub("DKI.DKI.","DKI.",dmri$feature)
dmri$feature=gsub("DTI.DTI.","DTI.",dmri$feature)
dmri$feature=gsub("SMT.SMT.","SMT.",dmri$feature)
dmri$feature=gsub("SMT_mc.SMT_mc.","SMT_mc.",dmri$feature)
dmri$feature=gsub("WMTI.WMTI.","WMTI.",dmri$feature)
## T1 features
T1$feature=gsub("lh_","",T1$feature)
T1$feature=gsub("Left.","",T1$feature)


# add a column of directionality for beta values / slopes
T1$Slope <- "No age relation (p > 0.05)"
T1$Slope[T1$beta_out > 0 & T1$log10p > -log10(0.05)] <- "T1w pos. related with age"
T1$Slope[T1$beta_out < 0 & T1$log10p > -log10(0.05)] <- "T1w neg. related with age"
##### THRESHOLDING TOP 10 FEATURES
data_new1 = T1[order(T1$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = T1[order(T1$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
T1$defeature = ifelse(T1$feature %in% top10$feature == TRUE, T1$feature, NA)

# old thresholding was based on beta values 
#T1$defeature = ifelse(abs(T1$beta_out) > 0.075, T1$feature, NA)

# check features, if interested
#filter(T1, is.na(defeature)==F)
# calculate log10 of the adjusted p-values
# select only relevant vars
T1 = T1 %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)
T1$defeature = gsub("Right.Lateral.Ventricle", "Lateral ventricle volume",T1$defeature)
T1$defeature = gsub("Right.Putamen", "Putamen volume",T1$defeature)
T1$defeature = gsub("rh.rostralmiddlefrontal_thickness", "Rostro-mid. thickness",T1$defeature)
T1$defeature = gsub("Right.Cerebellum.Cortex", "Cerebellum volume",T1$defeature)
T1$defeature = gsub("rhCorticalWhiteMatterVol", "WM volume",T1$defeature)
T1$defeature = gsub("Right.Amygdala", "Amygdala volume",T1$defeature)
T1$defeature = gsub("Right.Hippocampus", "Hippocampus volume",T1$defeature)
T1$defeature = gsub("rh.insula_area", "Insula area",T1$defeature)
T1$defeature = gsub("rh.caudalanteriorcingulate_thickness", "Caud. ant. cingulate thickness",T1$defeature)
T1$defeature = gsub("rh_WhiteSurfArea_area", "WM surface area",T1$defeature)
T1$defeature = gsub("Right.Inf.Lat.Vent", "Inferior lateral ventricle volume",T1$defeature)
T1$defeature = gsub("Right.Pallidum", "Pallidum volume",T1$defeature)
T1$defeature = gsub("Right.Thalamus.Proper", "Thalamus volume",T1$defeature)
T1$defeature = gsub("Right.Accumbens.area", "Accumbens volume",T1$defeature)
T1$defeature = gsub("rhCortexVol", "Cortex volume",T1$defeature)

# add a column of directionality for beta values / slopes
dmri$Slope <- "No age relation (p > 0.05)"
dmri$Slope[dmri$beta_out > 0 & dmri$log10p > -log10(0.05)] <- "dMRI pos. related with age"
dmri$Slope[dmri$beta_out < 0 & dmri$log10p > -log10(0.05)] <- "dMRI neg. related with age"

# # OLD
# # use only the thresholded features for labels (|beta| > 0.2)
# dmri$defeature = ifelse(abs(dmri$beta_out) > 0.2, dmri$feature, NA)
data_new1 = dmri[order(dmri$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = dmri[order(dmri$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
dmri$defeature = ifelse(dmri$feature %in% top10$feature == TRUE, dmri$feature, NA)

# display the features
na.omit(dmri)
# change these features (which will be labels) into a presentable format
dmri$defeature = gsub("DTI.FA_ILF", "DTI - FA ILF",dmri$defeature)
dmri$defeature = gsub("DTI.FA_SLTF", "DTI - FA SLTF",dmri$defeature)
dmri$defeature = gsub("SMT.smt_long_SLTF", "SMT - long SLTF",dmri$defeature)
dmri$defeature = gsub("DTI.RD_Fornix_Striaterminalis", "DTI - RD Fornix-Str.Term.",dmri$defeature)
dmri$defeature = gsub("DKI.RK_Superiorfrontooccipitalfasciculus", "DKI - RK Sup.fron.occ.Fasc.",dmri$defeature)
dmri$defeature = gsub("DTI.RD_Superiorfrontooccipitalfasciculus", "DTI - RD Sup.fron.occ.Fasc.",dmri$defeature)
dmri$defeature = gsub("SMT.smt_long_Tapetum", "SMT - long Tapetum",dmri$defeature)
dmri$defeature = gsub("SMT_mc.smt_mc_intra_Tapetum", "SMTmc - intra Tapetum",dmri$defeature)
dmri$defeature = gsub("BRIA.v_intra_Tapetum", "BRIA - Vintra Tapetum",dmri$defeature)
dmri$defeature = gsub("SMT.smt_fa_Tapetum", "SMT - FA Tapetum",dmri$defeature)
dmri$defeature = gsub("DTI.MD_Fornix_Striaterminalis", "DTI - MD Fornix-Str.Term.",dmri$defeature)
dmri$defeature = gsub("BRIA.v_csf_Fornix_Striaterminalis", "BRIA - vCSF Fornix-Str.Term.",dmri$defeature)
dmri$defeature = gsub("SMT_mc.smt_mc_extratrans_Cerebralpeduncle", "SMTmc - extratrans Cereb. peduncle",dmri$defeature)
dmri$defeature = gsub("BRIA.v_extra_Cerebralpeduncle", "BRIA - Vextra Cereb. peduncle",dmri$defeature)
dmri$defeature = gsub("SMT.smt_trans_Cerebralpeduncle", "SMT - trans Cereb. peduncle",dmri$defeature)
dmri$defeature = gsub("SMT_mc.smt_mc_extratrans_Superiorfrontooccipitalfasciculus", "SMTmc - extratrans Sup.fron.occ.Fasc.",dmri$defeature)
dmri$defeature = gsub("WMTI.axEAD_Superiorfrontooccipitalfasciculus", "WMTI - axEAD Sup.fron.occ.Fasc.",dmri$defeature)
dmri$defeature = gsub("WMTI.awf_Superiorfrontooccipitalfasciculus", "WMTI - AWF Sup.fron.occ.Fasc.",dmri$defeature)
dmri$defeature = gsub("BRIA.micro_Ax_SLTF", "BRIA - microAX SLTF",dmri$defeature)
dmri$defeature = gsub("BRIA.micro_Rd_CG", "BRIA - microRD Cingulate",dmri$defeature)
dmri$defeature = gsub("BRIA.micro_Rd_Fornix_Striaterminalis", "BRIA - microRD Fornix-Str.Term.",dmri$defeature)


# select relevant columns
dmri = dmri %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)

# merge dataframes
dat = rbind(T1,dmri)
dat$log10p = as.numeric(dat$log10p)
###### PLOT (both sexes) ####
plot = ggplot(data=dat, aes(x=beta_out, y=log10p,col = Slope, label=defeature)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(data = subset(dat, beta_out < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.6))+
  geom_text_repel(data = subset(dat, beta_out > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.3,0.6))+
  scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  xlab("Slope")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.5,0.5)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave(plot = plot, filename = "Volcano_Plot_abs.pdf", width = 13, height = 12)


####### PLOT PREP (males) ####
# remove prefixes and suffixes (remnant)
## dMRI features
dmri_m$feature=str_sub(dmri_m$feature,1,-2)
dmri_m$feature=gsub("BRIA.BRIA.","BRIA.",dmri_m$feature)
dmri_m$feature=gsub("DKI.DKI.","DKI.",dmri_m$feature)
dmri_m$feature=gsub("DTI.DTI.","DTI.",dmri_m$feature)
dmri_m$feature=gsub("SMT.SMT.","SMT.",dmri_m$feature)
dmri_m$feature=gsub("SMT_mc.SMT_mc.","SMT_mc.",dmri_m$feature)
dmri_m$feature=gsub("WMTI.WMTI.","WMTI.",dmri_m$feature)
## T1 features
T1_m$feature=gsub("lh_","",T1_m$feature)
T1_m$feature=gsub("Left.","",T1_m$feature)


# add a column of directionality for beta values / slopes
T1_m$Slope <- "No age relation (p > 0.05)"
T1_m$Slope[T1_m$beta_out > 0 & T1_m$log10p > -log10(0.05)] <- "T1w pos. related with age"
T1_m$Slope[T1_m$beta_out < 0 & T1_m$log10p > -log10(0.05)] <- "T1w neg. related with age"
##### THRESHOLDING TOP 10 FEATURES
data_new1 = T1_m[order(T1_m$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = T1_m[order(T1_m$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
T1_m$defeature = ifelse(T1_m$feature %in% top10$feature == TRUE, T1_m$feature, NA)

# old thresholding was based on beta values 
#T1_m$defeature = ifelse(abs(T1_m$beta_out) > 0.075, T1_m$feature, NA)

# check features, if interested
#filter(T1_m, is.na(defeature)==F)
# calculate log10 of the adjusted p-values
# select only relevant vars
T1_m = T1_m %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)
T1_m$defeature = gsub("Right.Lateral.Ventricle", "Lateral ventricle volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Putamen", "Putamen volume",T1_m$defeature)
T1_m$defeature = gsub("rh.rostralmiddlefrontal_thickness", "Rostro-mid. thickness",T1_m$defeature)
T1_m$defeature = gsub("Right.Cerebellum.Cortex", "Cerebellum volume",T1_m$defeature)
T1_m$defeature = gsub("rhCorticalWhiteMatterVol", "WM volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Amygdala", "Amygdala volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Hippocampus", "Hippocampus volume",T1_m$defeature)
T1_m$defeature = gsub("rh.insula_area", "Insula area",T1_m$defeature)
T1_m$defeature = gsub("rh.caudalanteriorcingulate_thickness", "Caud. ant. cingulate thickness",T1_m$defeature)
T1_m$defeature = gsub("rh_WhiteSurfArea_area", "WM surface area",T1_m$defeature)
T1_m$defeature = gsub("Right.Inf.Lat.Vent", "Inferior lateral ventricle volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Pallidum", "Pallidum volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Thalamus.Proper", "Thalamus volume",T1_m$defeature)
T1_m$defeature = gsub("Right.Accumbens.area", "Accumbens volume",T1_m$defeature)
T1_m$defeature = gsub("rhCortexVol", "Cortex volume",T1_m$defeature)
T1_m$defeature = gsub("rh.rostralanteriorcingulate_thickness", "Rost. ant. cingulate thickness",T1_m$defeature)
T1_m$defeature = gsub("Right.Caudate", "'Accumbens'Caudate volume",T1_m$defeature)

# add a column of directionality for beta values / slopes
dmri_m$Slope <- "No age relation (p > 0.05)"
dmri_m$Slope[dmri_m$beta_out > 0 & dmri_m$log10p > -log10(0.05)] <- "dMRI pos. related with age"
dmri_m$Slope[dmri_m$beta_out < 0 & dmri_m$log10p > -log10(0.05)] <- "dMRI neg. related with age"

# # OLD
# # use only the thresholded features for labels (|beta| > 0.2)
# dmri_m$defeature = ifelse(abs(dmri_m$beta_out) > 0.2, dmri_m$feature, NA)
data_new1 = dmri_m[order(dmri_m$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = dmri_m[order(dmri_m$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
dmri_m$defeature = ifelse(dmri_m$feature %in% top10$feature == TRUE, dmri_m$feature, NA)

# display the features
na.omit(dmri_m)
# change these features (which will be labels) into a presentable format
dmri_m$defeature = gsub("DTI.FA_ILF", "DTI - FA ILF",dmri_m$defeature)
dmri_m$defeature = gsub("DTI.FA_SLTF", "DTI - FA SLTF",dmri_m$defeature)
dmri_m$defeature = gsub("SMT.smt_long_SLTF", "SMT - long SLTF",dmri_m$defeature)
dmri_m$defeature = gsub("DTI.RD_Fornix_Striaterminalis", "DTI - RD Fornix-Str.Term.",dmri_m$defeature)
dmri_m$defeature = gsub("DKI.RK_Superiorfrontooccipitalfasciculus", "DKI - RK Sup.fron.occ.Fasc.",dmri_m$defeature)
dmri_m$defeature = gsub("DTI.RD_Superiorfrontooccipitalfasciculus", "DTI - RD Sup.fron.occ.Fasc.",dmri_m$defeature)
dmri_m$defeature = gsub("SMT.smt_long_Tapetum", "SMT - long Tapetum",dmri_m$defeature)
dmri_m$defeature = gsub("SMT_mc.smt_mc_intra_Tapetum", "SMTmc - intra Tapetum",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.v_intra_Tapetum", "BRIA - Vintra Tapetum",dmri_m$defeature)
dmri_m$defeature = gsub("SMT.smt_fa_Tapetum", "SMT - FA Tapetum",dmri_m$defeature)
dmri_m$defeature = gsub("DTI.MD_Fornix_Striaterminalis", "DTI - MD Fornix-Str.Term.",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.v_csf_Fornix_Striaterminalis", "BRIA - vCSF Fornix-Str.Term.",dmri_m$defeature)
dmri_m$defeature = gsub("SMT_mc.smt_mc_extratrans_Cerebralpeduncle", "SMTmc - extratrans Cereb. peduncle",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.v_extra_Cerebralpeduncle", "BRIA - Vextra Cereb. peduncle",dmri_m$defeature)
dmri_m$defeature = gsub("SMT.smt_trans_Cerebralpeduncle", "SMT - trans Cereb. peduncle",dmri_m$defeature)
dmri_m$defeature = gsub("SMT_mc.smt_mc_extratrans_Superiorfrontooccipitalfasciculus", "SMTmc - extratrans Sup.fron.occ.Fasc.",dmri_m$defeature)
dmri_m$defeature = gsub("WMTI.axEAD_Superiorfrontooccipitalfasciculus", "WMTI - axEAD Sup.fron.occ.Fasc.",dmri_m$defeature)
dmri_m$defeature = gsub("WMTI.awf_Superiorfrontooccipitalfasciculus", "WMTI - AWF Sup.fron.occ.Fasc.",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.micro_Ax_SLTF", "BRIA - microAX SLTF",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.micro_Rd_CG", "BRIA - microRD Cingulate",dmri_m$defeature)
dmri_m$defeature = gsub("BRIA.micro_Rd_Fornix_Striaterminalis", "BRIA - microRD Fornix-Str.Term.",dmri_m$defeature)
dmri_m$defeature = gsub("SMT_mc.smt_mc_extratrans_UF", "SMTmc - extratrans Uncinate Fasciculus",dmri_m$defeature)


# select relevant columns
dmri_m = dmri_m %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)

# merge dataframes
dat = rbind(T1_m,dmri_m)
dat$log10p = as.numeric(dat$log10p)
###### PLOT (males) ####
plot_m = ggplot(data=dat, aes(x=beta_out, y=log10p,col = Slope, label=defeature)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(data = subset(dat, beta_out < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.6))+
  geom_text_repel(data = subset(dat, beta_out > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.3,0.6))+
  scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  xlab("Slope")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.5,0.5)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave(plot = plot_m, filename = "Volcano_Plot_abs_males.pdf", width = 13, height = 12)


####### PLOT PREP (females) ####
# remove prefixes and suffixes (remnant)
## dmri_f features
dmri_f$feature=str_sub(dmri_f$feature,1,-2)
dmri_f$feature=gsub("BRIA.BRIA.","BRIA.",dmri_f$feature)
dmri_f$feature=gsub("DKI.DKI.","DKI.",dmri_f$feature)
dmri_f$feature=gsub("DTI.DTI.","DTI.",dmri_f$feature)
dmri_f$feature=gsub("SMT.SMT.","SMT.",dmri_f$feature)
dmri_f$feature=gsub("SMT_mc.SMT_mc.","SMT_mc.",dmri_f$feature)
dmri_f$feature=gsub("WMTI.WMTI.","WMTI.",dmri_f$feature)
## T1 features
T1_f$feature=gsub("lh_","",T1_f$feature)
T1_f$feature=gsub("Left.","",T1_f$feature)


# add a column of directionality for beta values / slopes
T1_f$Slope <- "No age relation (p > 0.05)"
T1_f$Slope[T1_f$beta_out > 0 & T1_f$log10p > -log10(0.05)] <- "T1w pos. related with age"
T1_f$Slope[T1_f$beta_out < 0 & T1_f$log10p > -log10(0.05)] <- "T1w neg. related with age"
##### THRESHOLDING TOP 10 FEATURES
data_new1 = T1_f[order(T1_f$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = T1_f[order(T1_f$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
T1_f$defeature = ifelse(T1_f$feature %in% top10$feature == TRUE, T1_f$feature, NA)

# old thresholding was based on beta values 
#T1_f$defeature = ifelse(abs(T1_f$beta_out) > 0.075, T1_f$feature, NA)

# check features, if interested
#filter(T1_f, is.na(defeature)==F)
# calculate log10 of the adjusted p-values
# select only relevant vars
T1_f = T1_f %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)
T1_f$defeature = gsub("Right.Lateral.Ventricle", "Lateral ventricle volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Putamen", "Putamen volume",T1_f$defeature)
T1_f$defeature = gsub("rh.rostralmiddlefrontal_thickness", "Rostro-mid. thickness",T1_f$defeature)
T1_f$defeature = gsub("Right.Cerebellum.Cortex", "Cerebellum volume",T1_f$defeature)
T1_f$defeature = gsub("rhCorticalWhiteMatterVol", "WM volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Amygdala", "Amygdala volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Hippocampus", "Hippocampus volume",T1_f$defeature)
T1_f$defeature = gsub("rh.insula_area", "Insula area",T1_f$defeature)
T1_f$defeature = gsub("rh.caudalanteriorcingulate_thickness", "Caud. ant. cingulate thickness",T1_f$defeature)
T1_f$defeature = gsub("rh_WhiteSurfArea_area", "WM surface area",T1_f$defeature)
T1_f$defeature = gsub("Right.Inf.Lat.Vent", "Inferior lateral ventricle volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Pallidum", "Pallidum volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Thalamus.Proper", "Thalamus volume",T1_f$defeature)
T1_f$defeature = gsub("Right.Accumbens.area", "Accumbens volume",T1_f$defeature)
T1_f$defeature = gsub("rhCortexVol", "Cortex volume",T1_f$defeature)

# add a column of directionality for beta values / slopes
dmri_f$Slope <- "No age relation (p > 0.05)"
dmri_f$Slope[dmri_f$beta_out > 0 & dmri_f$log10p > -log10(0.05)] <- "dmri_f pos. related with age"
dmri_f$Slope[dmri_f$beta_out < 0 & dmri_f$log10p > -log10(0.05)] <- "dmri_f neg. related with age"

# # OLD
# # use only the thresholded features for labels (|beta| > 0.2)
# dmri_f$defeature = ifelse(abs(dmri_f$beta_out) > 0.2, dmri_f$feature, NA)
data_new1 = dmri_f[order(dmri_f$beta_out, decreasing = T),]
largest = data_new1[1:5,]
data_new1 = dmri_f[order(dmri_f$beta_out, decreasing = F),]
lowest = data_new1[1:5,]
top10 = rbind(lowest, largest)
dmri_f$defeature = ifelse(dmri_f$feature %in% top10$feature == TRUE, dmri_f$feature, NA)

# display the features
na.omit(dmri_f)
# change these features (which will be labels) into a presentable format
dmri_f$defeature = gsub("DTI.FA_ILF", "DTI - FA ILF",dmri_f$defeature)
dmri_f$defeature = gsub("DTI.FA_SLTF", "DTI - FA SLTF",dmri_f$defeature)
dmri_f$defeature = gsub("SMT.smt_long_SLTF", "SMT - long SLTF",dmri_f$defeature)
dmri_f$defeature = gsub("DTI.RD_Fornix_Striaterminalis", "DTI - RD Fornix-Str.Term.",dmri_f$defeature)
dmri_f$defeature = gsub("DKI.RK_Superiorfrontooccipitalfasciculus", "DKI - RK Sup.fron.occ.Fasc.",dmri_f$defeature)
dmri_f$defeature = gsub("DTI.RD_Superiorfrontooccipitalfasciculus", "DTI - RD Sup.fron.occ.Fasc.",dmri_f$defeature)
dmri_f$defeature = gsub("SMT.smt_long_Tapetum", "SMT - long Tapetum",dmri_f$defeature)
dmri_f$defeature = gsub("SMT_mc.smt_mc_intra_Tapetum", "SMTmc - intra Tapetum",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.v_intra_Tapetum", "BRIA - Vintra Tapetum",dmri_f$defeature)
dmri_f$defeature = gsub("SMT.smt_fa_Tapetum", "SMT - FA Tapetum",dmri_f$defeature)
dmri_f$defeature = gsub("DTI.MD_Fornix_Striaterminalis", "DTI - MD Fornix-Str.Term.",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.v_csf_Fornix_Striaterminalis", "BRIA - vCSF Fornix-Str.Term.",dmri_f$defeature)
dmri_f$defeature = gsub("SMT_mc.smt_mc_extratrans_Cerebralpeduncle", "SMTmc - extratrans Cereb. peduncle",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.v_extra_Cerebralpeduncle", "BRIA - Vextra Cereb. peduncle",dmri_f$defeature)
dmri_f$defeature = gsub("SMT.smt_trans_Cerebralpeduncle", "SMT - trans Cereb. peduncle",dmri_f$defeature)
dmri_f$defeature = gsub("SMT_mc.smt_mc_extratrans_Superiorfrontooccipitalfasciculus", "SMTmc - extratrans Sup.fron.occ.Fasc.",dmri_f$defeature)
dmri_f$defeature = gsub("WMTI.axEAD_Superiorfrontooccipitalfasciculus", "WMTI - axEAD Sup.fron.occ.Fasc.",dmri_f$defeature)
dmri_f$defeature = gsub("WMTI.awf_Superiorfrontooccipitalfasciculus", "WMTI - AWF Sup.fron.occ.Fasc.",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.micro_Ax_SLTF", "BRIA - microAX SLTF",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.micro_Rd_CG", "BRIA - microRD Cingulate",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.micro_Rd_Fornix_Striaterminalis", "BRIA - microRD Fornix-Str.Term.",dmri_f$defeature)
dmri_f$defeature = gsub("BRIA.Dax_extra_SLTF", "BRIA - DAX extra SLTF",dmri_f$defeature)


# select relevant columns
dmri_f = dmri_f %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)

# merge dataframes
dat = rbind(T1_f,dmri_f)
dat$log10p = as.numeric(dat$log10p)
###### PLOT (females) ####
plot_f = ggplot(data=dat, aes(x=beta_out, y=log10p,col = Slope, label=defeature)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(data = subset(dat, beta_out < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.6))+
  geom_text_repel(data = subset(dat, beta_out > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.3,0.6))+
  scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  xlab("Slope")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.5,0.5)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave(plot = plot_f, filename = "Volcano_Plot_abs_females.pdf", width = 13, height = 12)