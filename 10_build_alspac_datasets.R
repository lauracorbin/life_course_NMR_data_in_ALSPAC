# code to build alspac format datasets

## input = metabolite data on aln/qlet IDs and processed feature sumstats and labels from 06
## output = data in alspac daatset format ready for pasisng to alspac data team for building into catalogue 

####################
###### SET UP ######
####################

# run code interactively from script directory

# capture output
rm(list=ls())
date <- Sys.Date()
sink(file=paste0("../../../results/10_build_datsets_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('dplyr')
library('ggplot2')
library('rlist')
library('data.table')
library('purrr')

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in linked data
data_list_new <- readRDS("../../../data/intermediate/07_linked_dataset_list.rds")
summary(data_list_new)

# read in feature metadata
feature_metadata <- readRDS("../../../data/intermediate/06_feature_sumstats_list.rds")
summary(feature_metadata)

# read in var labels
varlabels <- readRDS("../../../data/intermediate/06_var_labels_list.rds")
summary(varlabels)

# read in corstats results
corstats <- readRDS("../../../data/intermediate/08_corstats_list.rds")
summary(corstats)

# read in exclusions list (duplicates)
exclusions <- readRDS("../../../data/intermediate/09_exclusion_list.rds")
summary(exclusions)

####################################################################################
####################################################################################
## build G1 cohort dataset 
####################################################################################
####################################################################################

g1_data <- data_list_new[c(9,10,3,7,8,1,2)]
summary(g1_data)
g1_exclusions <- exclusions[c(9,10,3,7,8,1,2)]
summary(g1_exclusions)

# iterate over data from different timepoints
for (i in 1:7) {
  print(paste("processing data from:", names(g1_data[i])))
  # remove duplicates as identified in 09
  print("pre-exclude N:")
  print(nrow(g1_data[[i]]))
  g1_data[[i]] <- g1_data[[i]][!g1_data[[i]]$sample_id %in% g1_exclusions[[i]],]
  print("post-exclude N:")
  print(nrow(g1_data[[i]]))
  g1_data[[i]] <- as.data.frame(g1_data[[i]]) %>%
    select(aln, qlet, starts_with("has_sample"):ncol(.)) %>% 
    relocate(aln, qlet)
  g1_data[[i]]$aln <- as.integer(g1_data[[i]]$aln)
  g1_data[[i]]$qlet <- as.character(g1_data[[i]]$qlet)
  if (names(g1_data[i]) == "BBS_baseline") {
   names(g1_data)[i] <- "BBS_Baseline"
   names(g1_data[[i]])[3:ncol(g1_data[[i]])] <- paste0(names(g1_data[[i]])[3:ncol(g1_data[[i]])],"_Baseline")
  }
  if (names(g1_data[i]) == "BBS_postglu") {
    names(g1_data)[i] <- "BBS_PostGluc"
    names(g1_data[[i]])[3:ncol(g1_data[[i]])] <- paste0(names(g1_data[[i]])[3:ncol(g1_data[[i]])],"_PostGluc")
  }
}
summary(g1_data)

# Join datasets
g1_data_build <- g1_data %>% 
  reduce(full_join, by = c("aln", "qlet"))

# make those with no data to "no" for has sample (instead of NA)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_F7` = replace(`has_sample_at_F7`, is.na(has_sample_at_F7) , "no")) 
print("Samples at F7")
table(g1_data_build$has_sample_at_F7)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_F24` = replace(`has_sample_at_F24`, is.na(has_sample_at_F24) , "no"))
print("Samples at F24")
table(g1_data_build$has_sample_at_F24)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_F30` = replace(`has_sample_at_F30`, is.na(has_sample_at_F30) , "no")) 
print("Samples at F30")
table(g1_data_build$has_sample_at_F30)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_TF3` = replace(`has_sample_at_TF3`, is.na(has_sample_at_TF3) , "no")) 
print("Samples at TF3")
table(g1_data_build$has_sample_at_TF3)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_TF4` = replace(`has_sample_at_TF4`, is.na(has_sample_at_TF4) , "no")) 
print("Samples at TF4")
table(g1_data_build$has_sample_at_TF4)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_BBS_Baseline` = replace(`has_sample_at_BBS_Baseline`, is.na(has_sample_at_BBS_Baseline) , "no")) 
print("Samples at BBS Baseline")
table(g1_data_build$has_sample_at_BBS_Baseline)
g1_data_build <- g1_data_build %>% mutate(`has_sample_at_BBS_PostGluc` = replace(`has_sample_at_BBS_PostGluc`, is.na(has_sample_at_BBS_PostGluc) , "no")) 
print("Samples at BBS Baseline")
table(g1_data_build$has_sample_at_BBS_PostGluc)

####################################################################################
####################################################################################
## build G0 cohort dataset 
####################################################################################
####################################################################################

g0_data <- data_list_new[c(4,5)]
summary(g0_data)
g0_exclusions <- exclusions[c(4,5)]
summary(g0_exclusions)

# iterate over data from different timepoints
for (i in 1:2) {
  print(paste("processing data from:", names(g0_data[i])))
  # remove duplicates as identified in 09
  print("pre-exclude N:")
  print(nrow(g0_data[[i]]))
  g0_data[[i]] <- g0_data[[i]][!g0_data[[i]]$sample_id %in% g0_exclusions[[i]],]
  print("post-exclude N:")
  print(nrow(g0_data[[i]]))
  g0_data[[i]] <- as.data.frame(g0_data[[i]]) %>%
    select(aln, starts_with("has_sample"):ncol(.)) %>% 
    relocate(aln)
  g0_data[[i]]$aln <- as.integer(g0_data[[i]]$aln)
}
summary(g0_data)

# Join datasets
g0_data_build <- g0_data %>% 
  reduce(full_join, by = c("aln"))

# make those with no data to "no" for has sample (instead of NA)
print("Samples at FOM1")
g0_data_build <- g0_data_build %>% mutate(`has_sample_at_FOM1` = replace(`has_sample_at_FOM1`, is.na(has_sample_at_FOM1) , "no"))
table(g0_data_build$has_sample_at_FOM1)
print("Samples at FOM2")
g0_data_build <- g0_data_build %>% mutate(`has_sample_at_FOM2` = replace(`has_sample_at_FOM2`, is.na(has_sample_at_FOM2) , "no"))
table(g0_data_build$has_sample_at_FOM2)

####################################################################################
####################################################################################
## update variable labels - G1
####################################################################################
####################################################################################

g1_varlabels_orig <- c(varlabels[c(1)], list(NULL), varlabels[c(2,8,9,3,4)])
g1_varlabels <- c(varlabels[c(1)], list(NULL), varlabels[c(2,8,9,3,4)])
summary(g1_varlabels)
summary(g1_data)

for (i in 1:7){
  print(paste("processing data from:", names(g1_data[i])))
  # update BBS labels to include timepoint
  if (names(g1_data[i]) == "BBS_Baseline") {
    names(g1_varlabels)[i] = "BBS_Baseline"
    g1_varlabels[[i]] <- g1_varlabels[[i]] %>% 
      add_row(varname = "aln", varlabel = "aln", .before = 1)
    g1_varlabels[[i]][2,1] = "qlet"
    g1_varlabels[[i]][2,2] = "qlet"
    g1_varlabels[[i]]$varname[3:nrow(g1_varlabels[[i]])] <- paste0(g1_varlabels[[i]]$varname[3:nrow(g1_varlabels[[i]])], "_Baseline")
    g1_varlabels[[i]]$varlabel[3:nrow(g1_varlabels[[i]])] <- paste0(g1_varlabels[[i]]$varlabel[3:nrow(g1_varlabels[[i]])], "_Baseline")
    print(head(g1_varlabels[[i]]))
  }
  else if ((names(g1_data[i]) == "BBS_PostGluc") ) {
    names(g1_varlabels)[i] = "BBS_PostGluc"
    g1_varlabels[[i]] <- g1_varlabels_orig[[i-1]]
    g1_varlabels[[i]] <- g1_varlabels[[i]] %>% 
      add_row(varname = "aln", varlabel = "aln", .before = 1)
    g1_varlabels[[i]][2,1] = "qlet"
    g1_varlabels[[i]][2,2] = "qlet"
    g1_varlabels[[i]]$varname[3:nrow(g1_varlabels[[i]])] <- paste0(g1_varlabels[[i]]$varname[3:nrow(g1_varlabels[[i]])], "_PostGluc")
    g1_varlabels[[i]]$varlabel[3:nrow(g1_varlabels[[i]])] <- paste0(g1_varlabels[[i]]$varlabel[3:nrow(g1_varlabels[[i]])], "_PostGluc")
    print(head(g1_varlabels[[i]]))
  }
  else {
  g1_varlabels[[i]] <- g1_varlabels[[i]] %>% 
    add_row(varname = "aln", varlabel = "aln", .before = 1)
  g1_varlabels[[i]][2,1] = "qlet"
  g1_varlabels[[i]][2,2] = "qlet"
  print(head(g1_varlabels[[i]]))
  }
  # check all vars have labels 
  print(identical(g1_varlabels[[i]]$varname, names(g1_data[[i]])))
}

# add in row for dups indicator var (built in 07) at end of each varlabel file
g1_varlabels[["BBS_Baseline"]] <- g1_varlabels[["BBS_Baseline"]] %>% 
  add_row(varname = "had_dups_BBS_Baseline", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["BBS_PostGluc"]] <- g1_varlabels[["BBS_PostGluc"]] %>% 
  add_row(varname = "had_dups_BBS_PostGluc", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["F7"]] <- g1_varlabels[["F7"]] %>% 
  add_row(varname = "had_dups_F7", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["F24"]] <- g1_varlabels[["F24"]] %>% 
  add_row(varname = "had_dups_F24", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["F30"]] <- g1_varlabels[["F30"]] %>% 
  add_row(varname = "had_dups_F30", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["TF3"]] <- g1_varlabels[["TF3"]] %>% 
  add_row(varname = "had_dups_TF3", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)
g1_varlabels[["TF4"]] <- g1_varlabels[["TF4"]] %>% 
  add_row(varname = "had_dups_TF4", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=6)


# combine
g1_varlabels_all <- bind_rows(g1_varlabels)
dim(g1_varlabels_all)
head(g1_varlabels_all)
# remove duplicate aln/qlet rows
g1_varlabels_all <- unique(g1_varlabels_all)
dim(g1_varlabels_all)


## replace omiprep exclude labels that are too long for stata
g1_varlabels_all$varname <- gsub("omiprep_reason_exclude", "exclude_reason", g1_varlabels_all$varname)

####################################################################################
####################################################################################
## update variable labels - G0
####################################################################################
####################################################################################

g0_varlabels <- c(varlabels[c(5,6)])
summary(g0_varlabels)
summary(g0_data)

for (i in 1:2){
  print(paste("processing data from:", names(g0_data[i])))
  g0_varlabels[[i]][1,1] = "aln"
  g0_varlabels[[i]][1,2] = "aln"
  print(head(g0_varlabels[[i]]))
  # check all vars have labels 
  print(identical(g0_varlabels[[i]]$varname, names(g0_data[[i]])))
}


# add in row for dups indicator var (built in 07) at end of each varlabel file
g0_varlabels[["FOM1"]] <- g0_varlabels[["FOM1"]] %>% 
  add_row(varname = "had_dups_FOM1", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=5)
g0_varlabels[["FOM2"]] <- g0_varlabels[["FOM2"]] %>% 
  add_row(varname = "had_dups_FOM2", varlabel = "flag identifying participants who had a sample analysed from a re-invite", .after=5)

# combine
g0_varlabels_all <- bind_rows(g0_varlabels)
dim(g0_varlabels_all)
head(g0_varlabels_all)
# remove duplicate aln rows
g0_varlabels_all <- unique(g0_varlabels_all)
dim(g0_varlabels_all)

## replace omiprep exclude labels that are too long for stata
g0_varlabels_all$varname <- gsub("omiprep_reason_exclude", "exclude_reason", g0_varlabels_all$varname)


####################################################################################
####################################################################################
## write out datasets and var labels to return ALSPAC
####################################################################################
####################################################################################

# data 
#write.table(g1_data_build, paste0("../../../data/alspac_release/alnID/g1_updated_nmr_", date, ".csv"), sep=",", quote=F, row.names = F)
#write.table(g0_data_build, paste0("../../../data/alspac_release/alnID/g0_updated_nmr_", date, ".csv"), sep=",", quote=F, row.names = F)

# labels
write.table(g1_varlabels_all ,paste0("../../../data/alspac_release/alnID/g1_updated_nmr_labels_", date, ".txt"), sep="\t", quote=F, row.names = F)
write.table(g0_varlabels_all, paste0("../../../data/alspac_release/alnID/g0_updated_nmr_labels_", date, ".txt"), sep="\t", quote=F, row.names = F)



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
