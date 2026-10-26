# This script takes the raw data files and runs omiprep

## input = raw files
## output = metaboprep files
## run in RStudio from script folder

####################
###### SET UP ######
####################

# run code interactively from script directory
rm(list=ls())

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/05_omiprep_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('tibble')
library('omiprep')
library('kableExtra')
library('dendextend')
library('data.table')

####################################################################################
####################################################################################
## read in data to process
####################################################################################
####################################################################################

# feature metadata for all
feature_metadata <- read.table("../../../data/raw/metadata/feature_metadata.txt", sep="\t", h=T)
dim(feature_metadata)
# populate line for glycerol in feature metadata (not in annotation file used originally - from F30 - but does appear in BBS dataset)
glycerol <- c("Glycerol", "Glycerol", "Glycerol", "mmol/l", "Glycolysis related metabolites", NA, "Glycerol", "Glycerol (mmol/l)",
              "Glycolysis related metabolites", NA, "Glycerol (mmol/l)", "Glycerol", "no")
feature_metadata <- rbind(feature_metadata,glycerol)
dim(feature_metadata)
names(feature_metadata)[2] <- "feature_id"

##################
# bbs
bbs_data <- as.data.frame(fread("../../../data/raw/bbs/ALSPAC_ALBB-4_8_5.csv", sep=",",h=T))
dim(bbs_data)
bbs_data <- as.matrix(column_to_rownames(bbs_data, var = colnames(bbs_data)[1]))

bbs_samples <- read.table("../../../data/raw/bbs/ALSPAC_ALBB-4_8_5-sample-tags.csv", sep=",",h=T)
dim(bbs_samples)
names(bbs_samples)[1] <- "sample_id"

bbs <- Omiprep(data = bbs_data, samples = bbs_samples, features = feature_metadata)

##################
# f7
f7_data <- as.data.frame(fread("../../../data/raw/f7/ALSPAC_ALSPAC7-4_8_5.csv", sep=",",h=T))
dim(f7_data)
f7_data <- as.matrix(column_to_rownames(f7_data, var = colnames(f7_data)[1]))

f7_samples <- read.table("../../../data/raw/f7/ALSPAC_ALSPAC7-4_8_5-sample-tags.csv", sep=",",h=T)
dim(f7_samples)
names(f7_samples)[1] <- "sample_id"

f7 <- Omiprep(data = f7_data, samples = f7_samples, features = feature_metadata)

##################
# f24
f24_data <- as.data.frame(fread("../../../data/raw/f24/ALSPAC_ALTW-4_8_5.csv", sep=",",h=T))
dim(f24_data)
f24_data <- as.matrix(column_to_rownames(f24_data, var = colnames(f24_data)[1]))

f24_samples <- read.table("../../../data/raw/f24/ALSPAC_ALTW-4_8_5-sample-tags.csv", sep=",",h=T)
dim(f24_samples)
names(f24_samples)[1] <- "sample_id"

f24 <- Omiprep(data = f24_data, samples = f24_samples, features = feature_metadata)

##################
# f30
f30_data <- as.data.frame(fread("../../../data/raw/f30/f30-data.csv", sep=",",h=T))
dim(f30_data)
f30_data <- as.matrix(column_to_rownames(f30_data, var = colnames(f30_data)[1]))

f30_samples <- read.table("../../../data/raw/f30/f30-sample-tags.csv", sep=",",h=T)
dim(f30_samples)
names(f30_samples)[1] <- "sample_id"

f30 <- Omiprep(data = f30_data, samples = f30_samples, features = feature_metadata)

##################
# fom1
fom1_data <- as.data.frame(fread("../../../data/raw/fom1/ALSPAC_ALSPACM1-4_8_5.csv", sep=",",h=T))
dim(fom1_data)
fom1_data <- as.matrix(column_to_rownames(fom1_data, var = colnames(fom1_data)[1]))

fom1_samples <- read.table("../../../data/raw/fom1/ALSPAC_ALSPACM1-4_8_5-sample-tags.csv", sep=",",h=T)
dim(fom1_samples)
names(fom1_samples)[1] <- "sample_id"

fom1 <- Omiprep(data = fom1_data, samples = fom1_samples, features = feature_metadata)

##################
# fom2
fom2_data <- as.data.frame(fread("../../../data/raw/fom2/ALSPAC_ALSPACM2-4_8_5.csv", sep=",",h=T))
dim(fom2_data)
fom2_data <- as.matrix(column_to_rownames(fom2_data, var = colnames(fom2_data)[1]))

fom2_samples <- read.table("../../../data/raw/fom2/ALSPAC_ALSPACM2-4_8_5-sample-tags.csv", sep=",",h=T)
dim(fom2_samples)
names(fom2_samples)[1] <- "sample_id"

fom2 <- Omiprep(data = fom2_data, samples = fom2_samples, features = feature_metadata)

##################
# tf3
tf3_data <- as.data.frame(fread("../../../data/raw/tf3/ALSPAC_ALSPAC15-4_8_5.csv", sep=",",h=T))
dim(tf3_data)
tf3_data <- as.matrix(column_to_rownames(tf3_data, var = colnames(tf3_data)[1]))

tf3_samples <- read.table("../../../data/raw/tf3/ALSPAC_ALSPAC15-4_8_5-sample-tags.csv", sep=",",h=T)
dim(tf3_samples)
names(tf3_samples)[1] <- "sample_id"

tf3 <- Omiprep(data = tf3_data, samples = tf3_samples, features = feature_metadata)

##################
# tf4
tf4_data <- as.data.frame(fread("../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5_new_IDs.csv", sep=",",h=T))
dim(tf4_data)
tf4_data <- as.matrix(column_to_rownames(tf4_data, var = colnames(tf4_data)[1]))

tf4_samples <- read.table("../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-sample-tags_new_IDs.csv", sep=",",h=T)
dim(tf4_samples)
names(tf4_samples)[1] <- "sample_id"

tf4 <- Omiprep(data = tf4_data, samples = tf4_samples, features = feature_metadata)

##################
# g2
g2_data <- as.data.frame(fread("../../../data/raw/g2/ALSPAC_ALG2-4_8_5 id fix.tsv", sep="\t",h=T))
dim(g2_data)
g2_data <- as.matrix(column_to_rownames(g2_data, var = colnames(g2_data)[1]))

g2_samples <- read.table("../../../data/raw/g2/ALSPAC_ALG2-4_8_5-sample-tags_fixed.csv", sep=",",h=T)
dim(g2_samples)
names(g2_samples)[1] <- "sample_id"

g2 <- Omiprep(data = g2_data, samples = g2_samples, features = feature_metadata)

# dataset list
dataset_input_list <- list(bbs = bbs, f7 = f7, tf3 = tf3, tf4 = tf4, f24 = f24, f30 = f30, fom1 = fom1, fom2 = fom2, g2 = g2)

####################################################################################
####################################################################################
## run QC module
####################################################################################
####################################################################################

# make output list
dataset_output_list <- list(bbs = NA, f7 = NA, tf3 = NA, tf4 = NA, f24 = NA, f30 = NA, fom1 = NA, fom2 = NA, g2 = NA)

# run omiprep
for (i in 1:length(dataset_input_list)) {
  print(paste0("processing ", names(dataset_input_list)[i]))
  mydata <- dataset_input_list[[i]] |> quality_control(
    source_layer        = "input",
    sample_missingness  = 0.2,
    feature_missingness = 0.5,
    total_sum_abundance_sd  = NULL,
    outlier_udist       = 5,
    outlier_treatment   = "leave_be",
    pc_outlier_sd = NULL,
    tree_cut_height = 0.8
  )
  print(summary(mydata))
  dataset_output_list[[i]] <- mydata
}

# save out reports

generate_report(dataset_output_list[[1]], output_dir = "../../../data/raw/bbs/omiprep/", output_filename = "bbs")
generate_report(dataset_output_list[[2]], output_dir = "../../../data/raw/f7/omiprep/", output_filename = "f7", project = "f7")
generate_report(dataset_output_list[[3]], output_dir = "../../../data/raw/tf3/omiprep/", output_filename = "tf3", project = "tf3")
generate_report(dataset_output_list[[4]], output_dir = "../../../data/raw/tf4/omiprep/", output_filename = "tf4", project = "tf4")
generate_report(dataset_output_list[[5]], output_dir = "../../../data/raw/f24/omiprep/", output_filename = "f24", project = "f24")
generate_report(dataset_output_list[[6]], output_dir = "../../../data/raw/f30/omiprep/", output_filename = "f30", project = "f30")
generate_report(dataset_output_list[[7]], output_dir = "../../../data/raw/fom1/omiprep/", output_filename = "fom1", project = "fom1")
generate_report(dataset_output_list[[8]], output_dir = "../../../data/raw/fom2/omiprep/", output_filename = "fom2", project = "fom2")
generate_report(dataset_output_list[[9]], output_dir = "../../../data/raw/g2/omiprep/", output_filename = "g2", project = "g2")


# save out data
export(dataset_output_list[[1]], directory = "../../../data/raw/bbs/omiprep/", format = "omiprep")
export(dataset_output_list[[2]], directory = "../../../data/raw/f7/omiprep/", format = "omiprep")
export(dataset_output_list[[3]], directory = "../../../data/raw/tf3/omiprep/", format = "omiprep")
export(dataset_output_list[[4]], directory = "../../../data/raw/tf4/omiprep/", format = "omiprep")
export(dataset_output_list[[5]], directory = "../../../data/raw/f24/omiprep/", format = "omiprep")
export(dataset_output_list[[6]], directory = "../../../data/raw/f30/omiprep/", format = "omiprep")
export(dataset_output_list[[7]], directory = "../../../data/raw/fom1/omiprep/", format = "omiprep")
export(dataset_output_list[[8]], directory = "../../../data/raw/fom2/omiprep/", format = "omiprep")
export(dataset_output_list[[9]], directory = "../../../data/raw/g2/omiprep/", format = "omiprep")

# save as rds
saveRDS(dataset_output_list, "../../../data/intermediate/omiprep_data.rds")

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

