# Hemispheric brain age plot: beta coefficients by p-values
# Max Korbmacher, 02.06.2022, R version 4.1.2 (2021-11-01)

####### PREP ####

# load library
library(ggrepel)
library(ggpubr)
library(stringr)

# load data
T1 = read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_T1w_LI_pvals.csv")
dmri =  read.csv("/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/export7/Hemi_dMRI_LI_pvals.csv")

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
T1$Slope[T1$beta_out > 0 & T1$p_adj < 0.05] <- "T1w pos. related with age"
T1$Slope[T1$beta_out < 0 & T1$p_adj < 0.05] <- "T1w neg. related with age"
# use only thresholded (beta > 0.075) features for labels
T1$defeature = ifelse(abs(T1$beta_out) > 0.075, T1$feature, NA)
# check features, if interested
#filter(T1, is.na(defeature)==F)
# calculate log10 of the adjusted p-values
T1$log10p = -log10(T1$p_adj)
# select only relevant vars
T1 = T1 %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)
T1$defeature = gsub("Thalamus.Proper", "Thalamus volume",T1$defeature)
T1$defeature = gsub("Amygdala", "Amygdala volume",T1$defeature)
T1$defeature = gsub("Pallidum", "Pallidum volume",T1$defeature)
T1$defeature = gsub("Inf.Lat.Vent", "Inferior Lateral Ventricle volume",T1$defeature)
T1$defeature = gsub("_", " ",T1$defeature)
T1$defeature = gsub("inferiorparietal", "Inferior Parietal",T1$defeature)
T1$defeature = gsub("lateralorbitofrontal", "Lateral Orbitofrontal",T1$defeature)
T1$defeature = gsub("medialorbitofrontal", "Medial Orbitofrontal",T1$defeature)
T1$defeature = gsub("rostralmiddlefrontal", "Rostral Middle Frontal",T1$defeature)
T1$defeature = gsub("superiorfrontal", "Superior Frontal",T1$defeature)
T1$defeature = gsub("superiorparietal", "Superior Parietal",T1$defeature)
T1$defeature = gsub("entorhinal", "Entorhinal",T1$defeature)
T1$defeature = gsub("middletemporal", "Middle Temporal",T1$defeature)
T1$defeature = gsub("superiortemporal", "Superior Temporal",T1$defeature)
T1$defeature = gsub("WhiteSurfArea area", "WM Surface area",T1$defeature)

# add a column of directionality for beta values / slopes
dmri$Slope <- "No age relation (p > 0.05)"
dmri$Slope[dmri$beta_out > 0 & dmri$formatMpfr.logp. > -log10(0.05)] <- "dMRI pos. related with age"
dmri$Slope[dmri$beta_out < 0 & dmri$formatMpfr.logp. > -log10(0.05)] <- "dMRI neg. related with age"
# use only the thresholded features for labels (|beta| > 0.2)
dmri$defeature = ifelse(abs(dmri$beta_out) > 0.2, dmri$feature, NA)
# display the features
filter(dmri, is.na(defeature)==F)
# change these features (which will be labels) into a presentable format
dmri$defeature = gsub("SMT_mc.smt_mc_", "SMTmc - ",dmri$defeature)
dmri$defeature = gsub("BRIA.", "BRIA - ",dmri$defeature)
dmri$defeature = gsub("DTI.", "DTI - ",dmri$defeature)
dmri$defeature = gsub("SMT.smt_", "SMT - ",dmri$defeature)
dmri$defeature = gsub("Fornix_Striaterminalis", ": Fornix-Stria.",dmri$defeature)
dmri$defeature = gsub("_CG", ": Cingulum",dmri$defeature)
dmri$defeature = gsub("_UF", ": Uncin. Fasc.",dmri$defeature)
dmri$defeature = gsub("_Cerebralpeduncle", ": Cer. Ped.",dmri$defeature)
dmri$defeature = gsub("_", "",dmri$defeature)
# renmae log10p value vector
dmri$log10p = dmri$formatMpfr.logp.
# select relevant columns
dmri = dmri %>% dplyr::select(feature, beta_out, Slope, defeature, log10p)

# merge dataframes
dat = rbind(T1,dmri)

###### PLOT ####
plot = ggplot(data=dat, aes(x=beta_out, y=log10p,col = Slope, label=defeature)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(data = subset(dat, beta_out < -0.01),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.6))+
  geom_text_repel(data = subset(dat, beta_out > 0.01),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.3,0.6))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("Slope")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.5,0.5)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
plot


# 
# # plot T1
# T1plot = ggplot(data=T1, aes(x=beta_out, y=-log10(p_adj),col = Slope, label=defeature)) +
#   geom_point() + 
#   theme_minimal() +
#   geom_text_repel(
#     force        = 3.25,
#     #nudge_x      = -0.25,
#     direction    = "y",
#     hjust        = 1,
#     segment.size = 0.1,
#     max.overlaps = Inf)+
#   scale_color_manual(values=c("blue", "red", "black")) +
#   geom_vline(xintercept=c(-0.075, 0.075), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red") + xlab("Slope")+labs(title = "T1-weighted Features")
# 
# # plot dmri
# dmriplot = ggplot(data=dmri, aes(x=dmri$beta_out, y=-log10(p_adj),col = Slope, label=defeature)) +
#   geom_point() + 
#   theme_minimal() +
#   geom_hline(yintercept=-log10(min(dmri$p_adj)), linetype = "dashed",col="green")+
#   annotate("text", x = 0.3, y = max(-log10(dmri$p_adj)), color = "green", label = "Truncated y-axis \n appr. inf. (p~0)")+
#   #annotate("text", x = 0.3, y = max(-log10(dmri$p_adj))-5, color = "green", label = "")+
#   geom_text_repel(
#     force        = 3,
#     nudge_x      = -0.25,
#     direction    = "y",
#     hjust        = 1,
#     segment.size = 0.1,
#     max.overlaps = Inf
#   )+
#   scale_color_manual(values=c("blue", "red", "black")) +
#   #scale_color_viridis_d() +
#   geom_vline(xintercept=c(-0.2, 0.2), col="red") +
#   geom_hline(yintercept=-log10(0.05), col="red") + xlab("Slope") + labs(title="Diffusion Features")
# 
# plot = ggarrange(T1plot, dmriplot, nrow = 2, common.legend = T, legend = "bottom")
ggsave(plot = plot, filename = "/home/max/Documents/Projects/Brain_Age/Hemispheric_Brain_Age/Figures/Volcano_Plot.pdf", width = 12, height = 15)

