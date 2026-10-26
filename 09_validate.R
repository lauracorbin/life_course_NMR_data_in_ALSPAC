# code to validate against clinical vars and previous release

## input = metabolite data on aln/qlet IDs 
## output = correlation stats & figures, and list of duplicate sample IDs to exclude 

####################
###### SET UP ######
####################

# run code interactively from script directory

# capture output
rm(list=ls())
date <- Sys.Date()
sink(file=paste0("../../../results/09_validate_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('dplyr')
library('ggplot2')
library('blandr')
library('rlist')
library('data.table')
library('DescTools')

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in dataset list from 07
data_list_new <- readRDS("../../../data/intermediate/07_linked_dataset_list.rds")
summary(data_list_new)

# read in clinic measures
g1_measures <- read.table("../../../data/alspac_samples/g1_routine_measures.csv",h=T,sep=",")
g0_measures <- read.table("../../../data/alspac_samples/g0_routine_measures.csv",h=T,sep=",")

# read in original nmr release
g1_nmr <- as.data.frame(fread("../../../data/original_nmr_release/g1.csv",h=T,sep=","))
g0_nmr <- as.data.frame(fread("../../../data/original_nmr_release/g0m.csv",h=T,sep=","))

# read in multmum data
multmum <- read.csv("../../../data/multmum/multmum_data_clean.csv") %>%
  mutate(alnqlet = paste0(aln, "M"))

multmum_ids <- multmum %>%
  filter(!is.na(mz005m)) %>%
  .$alnqlet

# read in reinvite flag file
reinvites <- read.table("../../../data/linker_files/NMR_@30_F24_G2_re-invites.csv", sep=",", h=T)
dim(reinvites)
# check for dupslicate rows and delete
reinvites[duplicated(reinvites), ]
reinvites <- unique(reinvites)
dim(reinvites)

####################################################################################
####################################################################################
## reformat as needed
####################################################################################
####################################################################################

# recode missing values and keep only those with data available
# the last few rows in the g1_nmr file had a space in the aln_qlet column
# FOM1 NMR dataset has an additional -11 missing flag without label
g1_measures <- g1_measures %>%
  mutate(across(Glucose_BBS:Glucose_F30, as.numeric)) %>% # generates warning as missing labels converted to NA - this is intended
  filter(rowSums(!is.na(.[4:ncol(.)])) > 0)

g1_nmr <- g1_nmr %>%
  mutate(aln_qlet = paste0(aln, "_", qlet)) %>% 
  mutate(across(XXLVLDLP_F7:Gp_F24, as.numeric)) %>% # generates warning as missing labels converted to NA - this is intended
  filter(rowSums(!is.na(.[4:ncol(.)])) > 0)

g0_measures <- g0_measures %>%
  mutate(across(glucose_FOM1:trig_FOM2, as.numeric)) %>% # generates warning as missing labels converted to NA - this is intended
  filter(rowSums(!is.na(.[2:ncol(.)])) > 0)

g0_nmr <- g0_nmr %>%
  mutate(across(XXLVLDLP_FOM1:Sphingo_FOM2, as.numeric)) %>% # generates warning as missing labels converted to NA - this is intended
  mutate(across(XXLVLDLP_FOM1:Sphingo_FOM1, ~ ifelse(. == -11, NA, .))) %>%
  filter(rowSums(!is.na(.[3:ncol(.)])) > 0)

# add flags to clinical measures
names(g1_measures)[4:ncol(g1_measures)] <- paste0(names(g1_measures)[4:ncol(g1_measures)], "_clin")
names(g0_measures)[2:ncol(g0_measures)] <- paste0(names(g0_measures)[2:ncol(g0_measures)], "_clin")

####################################################################################
####################################################################################
## check that same individuals in new and old NMR releases
####################################################################################
####################################################################################

## f7
f7_ids_new <- unique(data_list_new$F7$alnqlet)
length(f7_ids_new)

f7_ids_old <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_F7),c("aln","qlet")]
f7_ids_old$alnqlet <- paste0(f7_ids_old$aln,f7_ids_old$qlet)
length(f7_ids_old$alnqlet)

length(which(!f7_ids_new %in% f7_ids_old$alnqlet))
length(which(!f7_ids_old$alnqlet %in% f7_ids_new))

print("IDs in updated dataset but not in original:")
f7_ids_new[which(!f7_ids_new %in% f7_ids_old$alnqlet)]
print("IDs in original dataset but not in updated:")
f7_ids_old$alnqlet[which(!f7_ids_old$alnqlet %in% f7_ids_new)]

#############################################################
## tf3
tf3_ids_new <- unique(data_list_new$TF3$alnqlet)
length(tf3_ids_new)
length(unique(tf3_ids_new))

tf3_ids_old <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_TF3),c("aln","qlet")]
tf3_ids_old$alnqlet <- paste0(tf3_ids_old$aln,tf3_ids_old$qlet)
length(tf3_ids_old$alnqlet)

length(which(!tf3_ids_new %in% tf3_ids_old$alnqlet))
length(which(!tf3_ids_old$alnqlet %in% tf3_ids_new))

print("IDs in updated dataset but not in original:")
tf3_ids_new[which(!tf3_ids_new %in% tf3_ids_old$alnqlet)]
print("IDs in original dataset but not in updated:")
tf3_ids_old$alnqlet[which(!tf3_ids_old$alnqlet %in% tf3_ids_new)]

#############################################################
## tf4
tf4_ids_new <- unique(data_list_new$TF4$alnqlet)
length(tf4_ids_new)
length(unique(tf4_ids_new))

tf4_ids_old <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_TF4),c("aln","qlet")]
tf4_ids_old$alnqlet <- paste0(tf4_ids_old$aln,tf4_ids_old$qlet)
length(tf4_ids_old$alnqlet)

length(which(!tf4_ids_new %in% tf4_ids_old$alnqlet))
length(which(!tf4_ids_old$alnqlet %in% tf4_ids_new))

print("IDs in updated dataset but not in original:")
tf4_ids_new[which(!tf4_ids_new %in% tf4_ids_old$alnqlet)]
print("IDs in original dataset but not in updated:")
tf4_ids_old$alnqlet[which(!tf4_ids_old$alnqlet %in% tf4_ids_new)]

#############################################################
## f24
f24_ids_new <- unique(data_list_new$F24$alnqlet)
length(f24_ids_new)
length(unique(f24_ids_new))

f24_ids_old <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_F24),c("aln","qlet")]
f24_ids_old$alnqlet <- paste0(f24_ids_old$aln,f24_ids_old$qlet)
length(f24_ids_old$alnqlet)

length(which(!f24_ids_new %in% f24_ids_old$alnqlet))
length(which(!f24_ids_old$alnqlet %in% f24_ids_new))

print("IDs in updated dataset but not in original:")
f24_ids_new[which(!f24_ids_new %in% f24_ids_old$alnqlet)]
print("IDs in original dataset but not in updated:")
f24_ids_old$alnqlet[which(!f24_ids_old$alnqlet %in% f24_ids_new)]


#############################################################
## f30
f30_ids_new <- unique(data_list_new$F30$alnqlet)
length(f30_ids_new)
length(unique(f30_ids_new))

#############################################################
## bbs
bbs_ids_new <- unique(c(data_list_new$BBS_baseline$alnqlet,data_list_new$BBS_postglu$alnqlet))
length(bbs_ids_new)
length(unique(bbs_ids_new))

bbs_ids_old_base <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_BBS_Baseline),c("aln","qlet")]
bbs_ids_old_end <- g1_nmr[!is.na(g1_nmr$XXLVLDLP_BBS_PostGluc),c("aln","qlet")]
bbs_ids_old <- rbind(bbs_ids_old_base,bbs_ids_old_end)
bbs_ids_old$alnqlet <- paste0(bbs_ids_old$aln,bbs_ids_old$qlet)
bbs_ids_old <- unique(bbs_ids_old$alnqlet)
length(bbs_ids_old)

length(which(!bbs_ids_new %in% bbs_ids_old))
length(which(!bbs_ids_old %in% bbs_ids_new))

print("IDs in updated dataset but not in original:")
bbs_ids_new[which(!bbs_ids_new %in% bbs_ids_old)]
print("IDs in original dataset but not in updated:")
bbs_ids_old[which(!bbs_ids_old %in% bbs_ids_new)]


#############################################################
## fom1
fom1_ids_new <- unique(data_list_new$FOM1$alnqlet)
length(fom1_ids_new)
length(unique(fom1_ids_new))

g0_nmr$alnqlet <- paste0(g0_nmr$aln,"M")
g0_nmr <- g0_nmr %>%
  relocate(alnqlet, .before = everything())
fom1_ids_old <- g0_nmr[!is.na(g0_nmr$XXLVLDLP_FOM1),c("alnqlet")]
#fom1_ids_old$alnqlet <- paste0(fom1_ids_old$aln,fom1_ids_old$qlet)
length(fom1_ids_old)

length(which(!fom1_ids_new %in% fom1_ids_old))
length(which(!fom1_ids_old %in% fom1_ids_new))

print("IDs in updated dataset but not in original:")
fom1_ids_new[which(!fom1_ids_new %in% fom1_ids_old)]
print("IDs in original dataset but not in updated:")
fom1_ids_old[which(!fom1_ids_old %in% fom1_ids_new)]
fom1_missing <- fom1_ids_old[which(!fom1_ids_old %in% fom1_ids_new)]
print("IDs missing from new FOM1 and not in mult mum list:")
fom1_missing[which(!fom1_missing %in% multmum_ids)]

#############################################################
## fom2
fom2_ids_new <- unique(data_list_new$FOM2$alnqlet)
length(fom2_ids_new)
length(unique(fom2_ids_new))

fom2_ids_old <- g0_nmr[!is.na(g0_nmr$XXLVLDLP_FOM2),c("alnqlet")]
length(fom2_ids_old)

length(which(!fom2_ids_new %in% fom2_ids_old))
length(which(!fom2_ids_old %in% fom2_ids_new))

print("IDs in updated dataset but not in original:")
fom2_ids_new[which(!fom2_ids_new %in% fom2_ids_old)]
print("IDs in original dataset but not in updated:")
fom2_ids_old[which(!fom2_ids_old %in% fom2_ids_new)]
fom2_missing <- fom2_ids_old[which(!fom2_ids_old %in% fom2_ids_new)]
print("IDs missing from new FOM1 and not in mult mum list:")
fom2_missing[which(!fom2_missing %in% multmum_ids)]

#######################################################################################
#######################################################################################
## check correlation/BA plot for clinical vs nmr 
#######################################################################################
#######################################################################################

#### G1 #### 

# process NMR data ahead of merge
g1_timepoints <- c("BBS_Baseline", "BBS_PostGluc", "F7", "F24", "TF3", "TF4", "F30")

g1_data <- data_list_new[c(9,10,3,1,7,8,2)]
summary(g1_data)
names(g1_data)[1] <- "BBS_Baseline"
names(g1_data)[2] <- "BBS_PostGluc"
summary(g1_data)

# processing for all timepoints
for (i in 1:7) {
  g1_data[[i]] <- g1_data[[i]] %>%
    mutate(aln_qlet = sub("^(.{5})", "\\1_", alnqlet)) %>% 
    relocate(aln_qlet, count_iids) %>% 
    select(aln_qlet, count_iids, starts_with("has_sample"):ncol(.))
  
  # Handle all renaming sequentially using base R names index
  names(g1_data[[i]])[2] <- paste0("count_iids_", names(g1_data)[i])
  if (names(g1_data)[i] == "BBS_Baseline"){
   names(g1_data[[i]])[3:ncol(g1_data[[i]])] <- paste0(names(g1_data[[i]])[3:ncol(g1_data[[i]])], "_Baseline")
  }
  if (names(g1_data)[i] == "BBS_PostGluc"){
    names(g1_data[[i]])[3:ncol(g1_data[[i]])] <- paste0(names(g1_data[[i]])[3:ncol(g1_data[[i]])], "_PostGluc")
  }
  # convert to dataframes
  g1_data[[i]] <- as.data.frame(g1_data[[i]])
}

# make list for with duplicates
g1_data_with_dups <- g1_data

# remove dups for validation
for (i in 1:7) {
  g1_data[[i]] <- g1_data_with_dups[[i]][g1_data_with_dups[[i]][,2] == 1,]
}
# check dims
lapply(g1_data, dim)
lapply(g1_data_with_dups, dim)

# merge g1 timepoints 
g1_nmr_new <- full_join(g1_data$F7, g1_data$BBS_Baseline, by = "aln_qlet") %>%
  full_join(., g1_data$BBS_PostGlu, by = "aln_qlet") %>%
  full_join(., g1_data$TF3, by = "aln_qlet") %>%
  full_join(., g1_data$TF4, by = "aln_qlet") %>%
  full_join(., g1_data$F24, by = "aln_qlet") %>%
  full_join(., g1_data$F30, by = "aln_qlet")

dim(g1_nmr_new)
length(unique(g1_nmr_new$aln_qlet))
length(unique(g1_nmr_new$aln_qlet))

# combine assayed glucose and nmr glucose
g1_glc_combined <- full_join(g1_measures, g1_nmr_new %>% select(aln_qlet, starts_with("Glucose_")), by = "aln_qlet")

# add assayed chol and nmr chol
g1_chol_combined <- full_join(g1_glc_combined, g1_nmr_new %>% select(aln_qlet, starts_with("Total_C_")), by = "aln_qlet")

# add assayed trig and nmr trig
g1_measures_combined <- full_join(g1_chol_combined, g1_nmr_new %>% select(aln_qlet, starts_with("Total_TG_")), by = "aln_qlet")

# reorder to do correlations
g1_measures_combined <- g1_measures_combined[,c("aln_qlet","Glucose_BBS_clin", "Glucose_BBS_Baseline",
                                             "pg_glucose_BBS_clin", "Glucose_BBS_PostGluc",
                                             "glucose_TF3_clin", "Glucose_TF3",
                                             "glucose_TF4_clin", "Glucose_TF4",
                                             "Glucose_F24_clin", "Glucose_F24",
                                             "Glucose_F30_clin", "Glucose_F30",
                                             "Chol_BBS_clin", "Total_C_BBS_Baseline", # no differentiate of pre/post
                                             "Chol_BBS_clin", "Total_C_BBS_PostGluc", # no differentiate of pre/post
                                             "CHOL_F7_clin",  "Total_C_F7",
                                             "chol_TF3_clin", "Total_C_TF3",
                                             "CHOL_TF4_clin", "Total_C_TF4",
                                             "Chol_F24_clin", "Total_C_F24",
                                             "Chol_F30_clin", "Total_C_F30",
                                             "Trig_BBS_clin", "Total_TG_BBS_Baseline", # no differentiate of pre/post
                                             "Trig_BBS_clin", "Total_TG_BBS_PostGluc", # no differentiate of pre/post
                                             "TRIG_F7_clin", "Total_TG_F7",
                                             "trig_TF3_clin", "Total_TG_TF3",
                                             "TRIG_TF4_clin", "Total_TG_TF4",
                                             "Trig_F24_clin", "Total_TG_F24")]

# calculate correlation between assayed measures and updated NMR 
plot_list = list()
plot_list_ba = list()

correlations <- data.frame(matrix(NA, nrow = 19, ncol = 6))
col_counts <- seq(2,ncol(g1_measures_combined),2)
count <- 0
plot_count <- 0
for (i in col_counts) {
  count = count+1
  plot_count = plot_count + 1
  correlations[count,1] <- names(g1_measures_combined)[i]
  correlations[count,2] <- names(g1_measures_combined)[i+1]
  correlations[count,3] <- sum(complete.cases(g1_measures_combined[,i:i+1]))
  correlations[count,4] <- SpearmanRho(g1_measures_combined[,i], g1_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[1]
  correlations[count,5] <- SpearmanRho(g1_measures_combined[,i], g1_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[2]
  correlations[count,6] <- SpearmanRho(g1_measures_combined[,i], g1_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[3]
  x_name <- names(g1_measures_combined)[i]
  y_name <- names(g1_measures_combined)[i+1]
  plot_list[[plot_count]] <- ggplot(g1_measures_combined, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
    labs(title = paste0("Updated NMR on clinical, r=", round(correlations[count,3],2)))
  plot_list_ba[[plot_count]] <- blandr.draw(g1_measures_combined[,i], g1_measures_combined[,i+1],
                                            plotTitle = paste0(names(g1_measures_combined)[i], "on ", names(g1_measures_combined)[i+1]) )
}
names(correlations) <- c("Clinic", "UpdatedNMR", "n", "Clinic_UpdatedNMR_corr", "LCI", "UCI")

print("G1: correlations between updated NMR and clinical/assayed measures")
print(correlations)

g1_correlations <- correlations

filename = "../../../results/figures/09_clinical_and_nmr_cor_G1.pdf"
outputfig <- ggpubr::ggarrange(plotlist=plot_list,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

filename = "../../../results/figures/09_clinical_and_nmr_ba_G1.pdf"
outputfig <- ggpubr::ggarrange(plotlist=plot_list_ba,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)


#### G0 #### 

# process NMR data ahead of merge
g0_timepoints <- c("FOM1", "FOM2")

g0_data <- data_list_new[c(4,5)]
summary(g0_data)

# processing for all timepoints
for (i in 1:2) {
  g0_data[[i]] <- g0_data[[i]] %>%
    relocate(aln, count_iids) %>% 
    select(aln, count_iids, starts_with("has_sample"):ncol(.))
  
  # Handle all renaming sequentially using base R names index
  names(g0_data[[i]])[2] <- paste0("count_iids_", names(g0_data)[i])
  # convert to dataframes
  g0_data[[i]] <- as.data.frame(g0_data[[i]])
}

# make list for with duplicates
g0_data_with_dups <- g0_data

# remove dups for validation
for (i in 1:2) {
  g0_data[[i]] <- g0_data_with_dups[[i]][g0_data_with_dups[[i]][,2] == 1,]
}
# check dims
lapply(g0_data, dim)
lapply(g0_data_with_dups, dim)

# merge g0 timepoints 
g0_nmr_new <- full_join(g0_data$FOM1, g0_data$FOM2, by = "aln") 

dim(g0_nmr_new)
length(unique(g0_nmr_new$aln))

# combine assayed glucose and nmr glucose
g0_nmr_new$aln <- as.numeric(g0_nmr_new$aln)
g0_glc_combined <- full_join(g0_measures, g0_nmr_new %>% select(aln, starts_with("Glucose_")), by = "aln")

# add assayed chol and nmr chol
g0_chol_combined <- full_join(g0_glc_combined, g0_nmr_new %>% select(aln, starts_with("Total_C_")), by = "aln")

# add assayed trig and nmr trig
g0_measures_combined <- full_join(g0_chol_combined, g0_nmr_new %>% select(aln, starts_with("Total_TG_")), by = "aln")

# reorder to do correlations
g0_measures_combined <- g0_measures_combined[,c("aln",
                                                "glucose_FOM1_clin", "Glucose_FOM1",
                                                "glucose_FOM2_clin", "Glucose_FOM2",
                                                "chol_FOM1_clin",  "Total_C_FOM1",
                                                "chol_FOM2_clin",  "Total_C_FOM2",
                                                "trig_FOM1_clin", "Total_TG_FOM1",
                                                "trig_FOM2_clin", "Total_TG_FOM2")]

# calculate correlation between assayed measures and updated NMR 
plot_list = list()
plot_list_ba = list()

correlations <- data.frame(matrix(NA, nrow = 6, ncol = 6))
col_counts <- seq(2,ncol(g0_measures_combined),2)
count <- 0
plot_count <- 0
for (i in col_counts) {
  count = count+1
  plot_count = plot_count + 1
  correlations[count,1] <- names(g0_measures_combined)[i]
  correlations[count,2] <- names(g0_measures_combined)[i+1]
  correlations[count,3] <- sum(complete.cases(g0_measures_combined[,i:i+1]))
  correlations[count,4] <- SpearmanRho(g0_measures_combined[,i], g0_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[1]
  correlations[count,5] <- SpearmanRho(g0_measures_combined[,i], g0_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[2]
  correlations[count,6] <- SpearmanRho(g0_measures_combined[,i], g0_measures_combined[,i+1], use = "pairwise.complete.obs", conf.level = 0.95)[3]
  x_name <- names(g0_measures_combined)[i]
  y_name <- names(g0_measures_combined)[i+1]
  plot_list[[plot_count]] <- ggplot(g0_measures_combined, aes(x = .data[[x_name]], y = .data[[y_name]])) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
    labs(title = paste0("Updated NMR on clinical, r=", round(correlations[count,3],2)))
  plot_list_ba[[plot_count]] <- blandr.draw(g0_measures_combined[,i], g0_measures_combined[,i+1],
                                            plotTitle = paste0(names(g0_measures_combined)[i], "on ", names(g0_measures_combined)[i+1]) )
}
names(correlations) <- c("Clinic", "UpdatedNMR", "n", "Clinic_UpdatedNMR_corr", "LCI", "UCI")

print("G0: correlations between updated NMR and clinical/assayed measures")
print(correlations)

g0_correlations <- correlations

filename = "../../../results/figures/09_clinical_and_nmr_cor_G0.pdf"
outputfig <- ggpubr::ggarrange(plotlist=plot_list,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)

filename = "../../../results/figures/09_clinical_and_nmr_ba_G0.pdf"
outputfig <- ggpubr::ggarrange(plotlist=plot_list_ba,
                               ncol=2, nrow=2,
                               common.legend = T,
                               legend="top",
                               heights = 2,
                               widths = 2)
ggpubr::ggexport(outputfig, filename=filename)


#######################################################################################
#######################################################################################
## check correlation of duplicates with clinical and retain one that matches the best
#######################################################################################

# G1
g1_data <- data_list_new[c(9,10,3,1,7,8,2)]
summary(g1_data)

g1_data_dups <- list(BBS_baseline = NULL, BBS_postglu = NULL, F7 = NULL, F24 = NULL, TF3 = NULL, TF4 = NULL, F30 = NULL)
g1_measures_combined_all <- list(BBS_baseline = NULL, BBS_postglu = NULL, F7 = NULL, F24 = NULL, TF3 = NULL, TF4 = NULL, F30 = NULL)

# merge clinical with nmr
for (i in 3:7) {
  print("processing:")
  print(names(g1_data[i]))
  g1_data[[i]] <- as.data.frame(g1_data[[i]])
  # identify dups alnqlet
  g1_data[[i]]$aln_qlet <- paste0(g1_data[[i]]$aln,"_",g1_data[[i]]$qlet) 
  dups <- g1_data[[i]]$alnqlet[which(g1_data[[i]]$count_iids == 2)]
  length(dups)
  g1_data_dups[[i]] <- g1_data[[i]][g1_data[[i]]$alnqlet %in% dups,]
  print(dim(g1_data_dups[[i]]))
  # merge in clinical data
  g1_measures_combined <- right_join(
    g1_measures, 
    g1_data_dups[[i]] %>% select(aln_qlet, count_iids, sample_id, alnqlet, starts_with(c("Glucose_", "Total_C_", "Total_TG_"))), 
    by = "aln_qlet", 
    relationship = "one-to-many"
  )
  # keep only those for this timepoint
  g1_measures_combined_all[[i]] <- g1_measures_combined %>% select(aln_qlet, count_iids, sample_id, alnqlet, contains(names(g1_data[i])))
  
}

## f7
data <- g1_measures_combined_all$F7
data$chol_diff <- abs(data$CHOL_F7_clin - data$Total_C_F7)
data$trig_diff <- abs(data$TRIG_F7_clin - data$Total_TG_F7)
data$total_diff <- data$chol_diff + data$trig_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln_qlet, .keep_all = TRUE)
# samples to drop
f7_dups_to_drop <-  data_unique$sample_id
length(f7_dups_to_drop)

## f24
data <- g1_measures_combined_all$F24
# merge on reinvite info 
data$reinvite <- reinvites$re.invites[match(data$sample_id,reinvites$sample_id)]
table(data$reinvite, useNA = "always")
# flag reinvites for removal
f24_dups_to_drop <- data$sample_id[which(data$reinvite == "True")]
length(f24_dups_to_drop)
# see if same samples are flagged with alternative approach
data$chol_diff <- abs(data$Chol_F24_clin - data$Total_C_F24)
data$trig_diff <- abs(data$Trig_F24_clin - data$Total_TG_F24)
data$glu_diff <- abs(data$Glucose_F24_clin - data$Glucose_F24)
data$total_diff <- data$chol_diff + data$trig_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln_qlet, .keep_all = TRUE)
# samples to drop
f24_dups_to_drop_test <-  data_unique$sample_id
length(f24_dups_to_drop_test)
length(which(f24_dups_to_drop %in% f24_dups_to_drop_test))

## f30
data <- g1_measures_combined_all$F30
# merge on reinvite info 
data$reinvite <- reinvites$re.invites[match(data$sample_id,reinvites$sample_id)]
table(data$reinvite, useNA = "always")
# flag reinvites for removal
f30_dups_to_drop <- data$sample_id[which(data$reinvite == "True")]
length(f30_dups_to_drop)
# see if same samples are flagged with alternative approach
data$chol_diff <- abs(data$Chol_F30_clin - data$Total_C_F30)
data$glu_diff <- abs(data$Glucose_F30_clin - data$Glucose_F30)
data$total_diff <- data$chol_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln_qlet, .keep_all = TRUE)
# samples to drop
f30_dups_to_drop_test <-  data_unique$sample_id
length(f30_dups_to_drop_test)
length(which(f30_dups_to_drop %in% f30_dups_to_drop_test))

## tf3
data <- g1_measures_combined_all$TF3
data$chol_diff <- abs(data$chol_TF3_clin - data$Total_C_TF3)
data$trig_diff <- abs(data$trig_TF3_clin - data$Total_TG_TF3)
data$glu_diff <- abs(data$glucose_TF3_clin - data$Glucose_TF3)
data$total_diff <- data$chol_diff + data$trig_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln_qlet, .keep_all = TRUE)
# samples to drop
tf3_dups_to_drop <-  data_unique$sample_id
length(tf3_dups_to_drop)


## tf4
data <- g1_measures_combined_all$TF4
data$chol_diff <- abs(data$CHOL_TF4_clin - data$Total_C_TF4)
data$trig_diff <- abs(data$TRIG_TF4_clin - data$Total_TG_TF4)
data$glu_diff <- abs(data$glucose_TF4_clin - data$Glucose_TF4)
data$total_diff <- data$chol_diff + data$trig_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln_qlet, .keep_all = TRUE)
# samples to drop
tf4_dups_to_drop <-  data_unique$sample_id
# add sample ID for single case where there were three samples analysed for an individual
tf4_dups_to_drop <- c(tf4_dups_to_drop, "**REDACTED**")
length(tf4_dups_to_drop)


# G0
g0_data <- data_list_new[c(4,5)]
summary(g0_data)

g0_data_dups <- list(FOM1 = NULL, FOM2 = NULL)
g0_measures_combined_all <- list(FOM1 = NULL, FOM2 = NULL)

# merge clinical with nmr
for (i in 1:2) {
  print("processing:")
  print(names(g0_data[i]))
  g0_data[[i]] <- as.data.frame(g0_data[[i]])
  # identify dups alnqlet
  #g0_data[[i]]$aln_qlet <- paste0(g0_data[[i]]$aln,"_",g0_data[[i]]$qlet) 
  dups <- g0_data[[i]]$aln[which(g0_data[[i]]$count_iids == 2)]
  length(dups)
  g0_data_dups[[i]] <- g0_data[[i]][g0_data[[i]]$aln %in% dups,]
  print(dim(g0_data_dups[[i]]))
  # merge in clinical data
  g0_data_dups[[i]]$aln <- as.integer(g0_data_dups[[i]]$aln)
  g0_measures_combined <- right_join(
    g0_measures, 
    g0_data_dups[[i]] %>% select(aln, count_iids, sample_id, starts_with(c("Glucose_", "Total_C_", "Total_TG_"))), 
    by = "aln", 
    relationship = "one-to-many"
  )
  # keep only those for this timepoint
  g0_measures_combined_all[[i]] <- g0_measures_combined %>% select(aln, count_iids, sample_id, contains(names(g0_data[i])))
  
}

## fom1
data <- g0_measures_combined_all$FOM1
data$chol_diff <- abs(data$chol_FOM1_clin - data$Total_C_FOM1)
data$trig_diff <- abs(data$trig_FOM1_clin - data$Total_TG_FOM1)
data$glu_diff <- abs(data$glucose_FOM1_clin - data$Glucose_FOM1)
data$total_diff <- data$chol_diff + data$trig_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln, .keep_all = TRUE)
# samples to drop
fom1_dups_to_drop <-  data_unique$sample_id
length(fom1_dups_to_drop)
# look at overlap with mult mums
data$alnqlet <- paste0(data$aln,"M")
data[which(data$alnqlet %in% multmum_ids),]

## fom2
data <- g0_measures_combined_all$FOM2
data$chol_diff <- abs(data$chol_FOM2_clin - data$Total_C_FOM2)
data$trig_diff <- abs(data$trig_FOM2_clin - data$Total_TG_FOM2)
data$glu_diff <- abs(data$glucose_FOM2_clin - data$Glucose_FOM2)
data$total_diff <- data$chol_diff + data$trig_diff + data$glu_diff
# sort with NAs at top
data <- data %>% 
  arrange(desc(is.na(total_diff)), desc(total_diff))
# id unique with biggest diff
data_unique <- data %>% 
  distinct(aln, .keep_all = TRUE)
# samples to drop
fom2_dups_to_drop <-  data_unique$sample_id
length(fom2_dups_to_drop)
# look at overlap with mult mums
data$alnqlet <- paste0(data$aln,"M")
data[which(data$alnqlet %in% multmum_ids),]

####################################################################################
####################################################################################
## save out 
####################################################################################
####################################################################################

# save out cor stats
write.table(g0_correlations, "../../../results/tables/09_G0_cor.csv", sep=",")
write.table(g1_correlations, "../../../results/tables/09_G1_cor.csv", sep=",")

# save out list of sample_ids to exclude
exclusion_list <- list(F24 = f24_dups_to_drop, F30 = f30_dups_to_drop, F7 = f7_dups_to_drop, FOM1 = fom1_dups_to_drop,
                      FOM2 = fom2_dups_to_drop, G2 = NULL, TF3 = tf3_dups_to_drop, TF4 = tf4_dups_to_drop, BBS_Baseline = NULL, BBS_PostGlu = NULL)
# save out
list.save(exclusion_list,"../../../data/intermediate/09_exclusion_list.rds")

##################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
rm(list=ls())
