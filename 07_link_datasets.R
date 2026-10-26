# This script takes as input the alspac ready files on lab IDs and links them to ALN/QLET provided by the ALSPAC data team

## input = processed metabolite data files & linker files
## output = metabolite data on aln/qlet IDs
## run in RStudio from script folder

####################
###### SET UP ######
####################

# capture output
rm(list=ls())
date <- Sys.Date()
sink(file=paste0("../../../results/07_link_datasets_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('dplyr')
library('rlist')
library('data.table')

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in data files
data_file_names <- list.files(path = "../../../data/alspac_release/labID", pattern="nmr_data", full.names=TRUE)

data_files_list <- list(BBS=NULL, F24 = NULL, F30 = NULL, F7 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(data_file_names)){
  data_files_list[[i]] <- as.data.frame(fread(file=data_file_names[[i]], h=T, sep="\t"))
}

# read in linker files
f30_linker <- read.table("../../../data/linker_files/NMRat30_idmatched_ALN.csv", sep=",", h=T)
print("dimensions of f30 linker")
dim(f30_linker)

main_linker <- read.table("../../../data/linker_files/idswap_from_LF/idswap_v4.csv", sep=",", h=T)
print("dimensions of main linker")
dim(main_linker)

####################################################################################
####################################################################################
## switch sample IDs for ALN/QLET - check linkages (unlinked and duplicated)
####################################################################################
####################################################################################

### create aln and qlet separate in linker file
main_linker$aln <- gsub("[a-zA-Z]","",main_linker$alnqlet)
main_linker$qlet <- gsub("[0-9]","",main_linker$alnqlet)

### F30
print("processing f30")
# merge linker with data - keep all rows of both (to check for missing IDs and/or missing data)
print("dimensions of metabolite data")
dim(data_files_list$F30)
print("dimensions of linker file")
dim(f30_linker)
f30_linker <- distinct(f30_linker) # for some reason, two lines are duplicated in the F30 linker file. Removed them here.
print("dimensions of linker file after excluding duplicated lines")
dim(f30_linker)
print("count of unique sample ids in linker file")
length(unique(f30_linker$sample_id))
print("count of unique aln/qlet ids in linker file")
length(unique(paste0(f30_linker$aln,f30_linker$qlet)))
f30_linker$alnqlet <- paste0(f30_linker$aln,f30_linker$qlet)
# merge linker and data
f30 <- full_join(f30_linker,data_files_list$F30, by = c("raw_ids" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(f30)

# count no. that appear in linker but don't have data
print("Check for samples in linker file but not in data file")
table(f30$has_sample_at, useNA="always")

# count samples not linked
print("No. of smaples that are not linked to aln/qlet")
length(which(is.na(f30$alnqlet)))

# check uniqueness of IDs
print("count of unique 'raw_ids' (an extension of sample_id)")
length(unique(f30$raw_ids))
print("count of unique sample ids")
length(unique(f30$sample_id))
print("count of unique aln/qlet ids")
length(unique(f30$alnqlet))

f30 <- f30 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(f30$count_iids)

# build aln, qlet dataset
data_files_list[["F30"]] <- f30 |> relocate(alnqlet, count_iids)

### F7
print("processing f7")
print("dimensions of metabolite data")
dim(data_files_list$F7)
f7_linker <- main_linker %>% filter(study_type == "ALSPAC7")
print("dimensions of linker file")
dim(f7_linker)
print("No. of unique lines in linker file")
n_distinct(f7_linker)
print("count of unique 'individual_ids' in linker file")
length(unique(f7_linker$individual_id))
print("count of unique aln/qlet ids in linker file")
length(unique(f7_linker$alnqlet))
print("count of unique sample ids in linker file")
length(unique(f7_linker$sample_id))
print("count of unique sample codes in data file")
length(unique(data_files_list$F7$sample_id))
# merge linker and data
f7 <- full_join(f7_linker, data_files_list$F7, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(f7)
# count no. that appear in linker but don't have data
print("Check for samples in linker file but not in data file")
table(f7$has_sample_at, useNA="always")
# drop those with no data
f7 <- f7[!is.na(f7$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(f7)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(f7$alnqlet)))
# count repeated iids
f7 <- f7 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(f7$count_iids)

# build aln, qlet dataset
data_files_list[["F7"]] <- f7 |> relocate(alnqlet, count_iids)

### TF3
print("processing tf3")
print("dimensions of metabolite data")
dim(data_files_list$TF3)
tf3_linker <- main_linker %>% filter(study_type == "ALSPAC15")
print("dimensions of linker file")
dim(tf3_linker)
print("No. of unique lines in linker file")
n_distinct(tf3_linker)
print("count of unique 'individual_ids'")
length(unique(tf3_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(tf3_linker$alnqlet))
print("count of unique sample ids")
length(unique(tf3_linker$sample_id))
print("count of unique sample codes")
length(unique(data_files_list$TF3$sample_id))
tf3 <- full_join(tf3_linker, data_files_list$TF3, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(tf3)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(tf3$has_sample_at, useNA="always")
# drop those with no data
tf3 <- tf3[!is.na(tf3$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(tf3)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(tf3$alnqlet)))
# count repeats
tf3 <- tf3 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(tf3$count_iids)

# build aln, qlet dataset
data_files_list[["TF3"]] <- tf3 |> relocate(alnqlet, count_iids)

### TF4
print("processing tf4")
print("dimensions of metabolite data")
dim(data_files_list$TF4)
tf4_linker <- main_linker %>% filter(study_type == "ALSPAC17")
print("dimensions of linker file")
dim(tf4_linker)
tf4_linker <- tf4_linker %>%
  mutate(sample_id = paste(sample_id, row_number(), sep="_"))
print("No. of unique lines in linker file")
n_distinct(tf4_linker)
print("count of unique 'individual_ids'")
length(unique(tf4_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(tf4_linker$alnqlet))
print("count of unique sample ids")
length(unique(tf4_linker$sample_id))
print("count of unique sample codes in linker")
length(unique(tf4_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$TF4$sample_id))
tf4 <- full_join(tf4_linker, data_files_list$TF4, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(tf4)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(tf4$has_sample_at, useNA="always")
# drop those with no data
tf4 <- tf4[!is.na(tf4$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(tf4)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(tf4$alnqlet)))
# count repeated iids
tf4 <- tf4 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(tf4$count_iids)

# build aln, qlet dataset
data_files_list[["TF4"]] <- tf4 |> relocate(alnqlet, count_iids)

### F24
print("processing f24")
print("dimensions of metabolite data")
dim(data_files_list$F24)
f24_linker <- main_linker %>% filter(study_type == "ALTW")
print("dimensions of linker file")
dim(f24_linker)
print("No. of unique lines in linker file")
n_distinct(f24_linker)
print("count of unique 'individual_ids'")
length(unique(f24_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(f24_linker$alnqlet))
print("count of unique sample ids")
length(unique(f24_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$F24$sample_id))
f24 <- full_join(f24_linker, data_files_list$F24, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(f24)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(f24$has_sample_at, useNA="always")
# drop those with no data
f24 <- f24[!is.na(f24$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(f24)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(f24$alnqlet)))
# count repeated iids
f24 <- f24 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(f24$count_iids)

# build aln, qlet dataset
data_files_list[["F24"]] <- f24 |> relocate(alnqlet, count_iids)


### FOM1
print("processing fom1")
print("dimensions of metabolite data")
dim(data_files_list$FOM1)
fom1_linker <- main_linker %>% filter(study_type == "ALSPACM1")
print("dimensions of linker file")
dim(fom1_linker)
print("No. of unique lines in linker file")
n_distinct(fom1_linker)
# find missing aln/qlet and remove
na_row <- which(fom1_linker$alnqlet == "nanM")
fom1_linker <- fom1_linker[-na_row,]
print("dimensions of linker file after NA row deleted")
dim(fom1_linker)
print("count of unique 'individual_ids'")
length(unique(fom1_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(fom1_linker$alnqlet))
print("count of unique sample ids")
length(unique(fom1_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$FOM1$sample_id))
fom1 <- full_join(fom1_linker, data_files_list$FOM1, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(fom1)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(fom1$has_sample_at, useNA="always")
# drop those with no data
fom1 <- fom1[!is.na(fom1$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(fom1)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(fom1$alnqlet)))
# drop unlinked
fom1 <- fom1[!is.na(fom1$alnqlet),]
print("dimensions of merged metabolite data with linker file after removing those with no data and those not linked")
dim(fom1)
# count repeated iids
fom1 <- fom1 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(fom1$count_iids) # 226 rows without iids but linked to alnqlet

# build aln, qlet dataset
data_files_list[["FOM1"]] <- fom1 |> relocate(alnqlet, count_iids)

### FOM2
print("processing fom2")
print("dimensions of metabolite data")
dim(data_files_list$FOM2)
fom2_linker <- main_linker %>% filter(study_type == "ALSPACM2")
print("dimensions of linker file")
dim(fom2_linker)
print("No. of unique lines in linker file")
n_distinct(fom2_linker)
print("count of unique 'individual_ids'")
length(unique(fom2_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(fom2_linker$alnqlet))
print("count of unique sample ids")
length(unique(fom2_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$FOM2$sample_id))
fom2 <- full_join(fom2_linker, data_files_list$FOM2, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(fom2)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(fom2$has_sample_at, useNA="always")
# drop those with no data
fom2 <- fom2[!is.na(fom2$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(fom2)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(fom2$alnqlet)))
# count repeated iids
fom2 <- fom2 %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(fom2$count_iids)

# build aln, qlet dataset
data_files_list[["FOM2"]] <- fom2 |> relocate(alnqlet, count_iids)

### G2
print("processing g2")
print("dimensions of metabolite data")
dim(data_files_list$G2)
g2_linker <- main_linker %>% filter(study_type == "ALG2")
print("dimensions of linker file")
dim(g2_linker)
print("No. of unique lines in linker file")
n_distinct(g2_linker)
print("count of unique 'individual_ids'")
length(unique(g2_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(g2_linker$alnqlet))
print("count of unique sample ids")
length(unique(g2_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$G2$sample_id))
g2 <- full_join(g2_linker, data_files_list$G2, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(g2)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(g2$has_sample_at, useNA="always")
# drop those with no data
g2 <- g2[!is.na(g2$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(g2)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(g2$alnqlet)))
# count repeated iids
g2 <- g2 %>%
  group_by(individual_id) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(g2$count_iids)

# build aln, qlet dataset
data_files_list[["G2"]] <- g2 |> relocate(alnqlet, count_iids)

### BBS
print("processing bbs")
print("dimensions of metabolite data")
dim(data_files_list$BBS)
bbs_linker <- main_linker %>% filter(study_type == "ALBB")
print("dimensions of linker file")
dim(bbs_linker)
print("No. of unique lines in linker file")
n_distinct(bbs_linker)
print("count of unique 'individual_ids'")
length(unique(bbs_linker$individual_id))
print("count of unique aln/qlet ids")
length(unique(bbs_linker$alnqlet))
print("count of unique sample ids")
length(unique(bbs_linker$sample_id))
print("count of unique sample codes in data")
length(unique(data_files_list$BBS$sample_id))
bbs <- full_join(bbs_linker, data_files_list$BBS, by = c("sample_id" = "sample_id"))
print("dimensions of merged metabolite data with linker file")
dim(bbs)
# count no. that appear in linker but dont have data
print("Check for samples in linker file but not in data file")
table(bbs$has_sample_at, useNA="always")
# drop those with no data
bbs <- bbs[!is.na(bbs$has_sample_at),]
print("dimensions of merged metabolite data with linker file after removing those with no data")
dim(bbs)
# count samples not linked
print("No. of samples that are not linked to aln/qlet")
length(which(is.na(bbs$alnqlet)))
# count repeated iids
bbs <- bbs %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id")
table(bbs$count_iids)

# Split BBS dataset into baseline and post-glucose for ease of renaming variables
bbs_baseline <- bbs %>%
  filter(grepl("415$", sample_id))
bbs_postgluc <- bbs %>%
  filter(grepl("420$", sample_id))

# count repeated iids
bbs_baseline <- bbs_baseline %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id (baseline)")
table(bbs_baseline$count_iids)

bbs_postgluc <- bbs_postgluc %>%
  group_by(alnqlet) %>%
  mutate(count_iids = row_number()) %>%
  ungroup() # Ungroup to prevent unexpected behavior downstream
print("Count no. of duplicates by aln/qlet id (postglu)")
table(bbs_postgluc$count_iids)

# build aln, qlet dataset
data_files_list[["BBS_baseline"]] <- bbs_baseline |> relocate(alnqlet, count_iids)
data_files_list[["BBS_postglu"]] <- bbs_postgluc |> relocate(alnqlet, count_iids)

# drop BBS combined
data_files_list[["BBS"]] <- NULL

####################################################################################
####################################################################################
## derive var to capture which individuals had repeat samples
####################################################################################
####################################################################################


## G1
f7_dups <- f7[f7$alnqlet %in% (f7[f7$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(f7_dups)
f7_nondups <- f7[!f7$alnqlet %in% (f7[f7$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(f7_dups) + nrow(f7_nondups) != nrow(f7)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
f7_dups$had_duplicate_samples_F7 <- "yes"
f7_nondups$had_duplicate_samples_F7 <- "no"
f7_dups_var <- rbind(f7_dups,f7_nondups)
table(f7_dups_var$had_duplicate_samples_F7)

f24_dups <- f24[f24$alnqlet %in% (f24[f24$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(f24_dups)
f24_nondups <- f24[!f24$alnqlet %in% (f24[f24$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(f24_dups) + nrow(f24_nondups) != nrow(f24)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
f24_dups$had_duplicate_samples_F24 <- "yes"
f24_nondups$had_duplicate_samples_F24 <- "no"
f24_dups_var <- rbind(f24_dups,f24_nondups)
table(f24_dups_var$had_duplicate_samples_F24)

f30_dups <- f30[f30$alnqlet %in% (f30[f30$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(f30_dups)
f30_nondups <- f30[!f30$alnqlet %in% (f30[f30$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(f30_dups) + nrow(f30_nondups) != nrow(f30)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
f30_dups$had_duplicate_samples_F30 <- "yes"
f30_nondups$had_duplicate_samples_F30 <- "no"
f30_dups_var <- rbind(f30_dups,f30_nondups)
table(f30_dups_var$had_duplicate_samples_F30)

tf3_dups <- tf3[tf3$alnqlet %in% (tf3[tf3$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(tf3_dups)
tf3_nondups <- tf3[!tf3$alnqlet %in% (tf3[tf3$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(tf3_dups) + nrow(tf3_nondups) != nrow(tf3)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
tf3_dups$had_duplicate_samples_TF3 <- "yes"
tf3_nondups$had_duplicate_samples_TF3 <- "no"
tf3_dups_var <- rbind(tf3_dups,tf3_nondups)
table(tf3_dups_var$had_duplicate_samples_TF3)

tf4_dups <- tf4[tf4$alnqlet %in% (tf4[tf4$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(tf4_dups)
tf4_nondups <- tf4[!tf4$alnqlet %in% (tf4[tf4$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(tf4_dups) + nrow(tf4_nondups) != nrow(tf4)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
tf4_dups$had_duplicate_samples_TF4 <- "yes"
tf4_nondups$had_duplicate_samples_TF4 <- "no"
tf4_dups_var <- rbind(tf4_dups,tf4_nondups)
table(tf4_dups_var$had_duplicate_samples_TF4)

bbs_baseline_dups_var <- bbs_baseline[!duplicated(bbs_baseline$alnqlet),c("alnqlet")]
dim(bbs_baseline_dups_var)
bbs_baseline_dups_var$had_duplicate_samples_BBS_Baseline <- "no"
table(bbs_baseline_dups_var$had_duplicate_samples_BBS_Baseline, useNA = "always")

bbs_postgluc_dups_var <- bbs_postgluc[!duplicated(bbs_postgluc$alnqlet),c("alnqlet")]
dim(bbs_postgluc_dups_var)
bbs_postgluc_dups_var$had_duplicate_samples_BBS_PostGluc <- "no"
table(bbs_postgluc_dups_var$had_duplicate_samples_BBS_PostGluc, useNA = "always")

## combine into one list
df_list <- list(f7_dups_var, tf3_dups_var, tf4_dups_var, f24_dups_var, f30_dups_var, bbs_baseline_dups_var, bbs_postgluc_dups_var)
g1_dups_vars <- df_list %>% reduce(full_join, by = "alnqlet")
dim(g1_dups_vars)


## G0

fom1_dups <- fom1[fom1$alnqlet %in% (fom1[fom1$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(fom1_dups)
fom1_nondups <- fom1[!fom1$alnqlet %in% (fom1[fom1$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(fom1_dups) + nrow(fom1_nondups) != nrow(fom1)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
fom1_dups$had_duplicate_samples_FOM1 <- "yes"
fom1_nondups$had_duplicate_samples_FOM1 <- "no"
fom1_dups_var <- rbind(fom1_dups,fom1_nondups)
table(fom1_dups_var$had_duplicate_samples_FOM1)

fom2_dups <- fom2[fom2$alnqlet %in% (fom2[fom2$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
dim(fom2_dups)
fom2_nondups <- fom2[!fom2$alnqlet %in% (fom2[fom2$count_iids != 1, c("alnqlet")])$alnqlet, c("alnqlet")]
# Standard if statement validation
if (nrow(fom2_dups) + nrow(fom2_nondups) != nrow(fom2)) {
  stop("Validation failed: val1 and val2 do not add up to total.")
}
fom2_dups$had_duplicate_samples_FOM2 <- "yes"
fom2_nondups$had_duplicate_samples_FOM2 <- "no"
fom2_dups_var <- rbind(fom2_dups,fom2_nondups)
table(fom2_dups_var$had_duplicate_samples_FOM2)

## combine into one list
df_list <- list(fom1_dups_var,fom2_dups_var)
g0_dups_vars <- df_list %>% reduce(full_join, by = "alnqlet")
dim(g0_dups_vars)

####################################################################################
####################################################################################
## save out 
####################################################################################
####################################################################################

# save out
#list.save(data_files_list,"../../../data/intermediate/07_linked_dataset_list.rds")

# save out list of duplicates
f7_dups <- f7[f7$alnqlet %in% (f7[f7$count_iids == 2, c("alnqlet")])$alnqlet, c(1:7, ncol(f7))]
f7_dups <- f7_dups %>% arrange(alnqlet, count_iids)
write.table(f7_dups, "../../../data/intermediate/07_f7dups.txt", row.names=F, quote=F)
f7_nondups <- f7[!f7$alnqlet %in% (f7[f7$count_iids == 2, c("alnqlet")])$alnqlet, c(1:7, ncol(f7))]
f7_nondups <- f7_nondups %>% arrange(alnqlet)
write.table(f7_nondups, "../../../data/intermediate/07_f7nondups.txt", row.names=F, quote=F)

# save out variable that flags where individuals had duplicate samples analysed (because of re-invites)
write.table(g0_dups_vars,"../../../data/alspac_release/alnID/g0_duplicate_samples_var.csv", sep=",",row.names=F,quote=F)
write.table(g1_dups_vars,"../../../data/alspac_release/alnID/g1_duplicate_samples_var.csv", sep=",",row.names=F,quote=F)

######################################################
## QUIT SESSION ##
##################

# capture session info
print("Session information:")
sessionInfo()

####################################################################################
# save output
sink()
rm(list=ls())
