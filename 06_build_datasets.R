# This script takes the output from metaboprep and builds data into a format suitable for ALSPAC

## input = omiprep output files
## output = alspac build to return for linking
## run in RStudio from script folder

####################
###### SET UP ######
####################

# run code interactively from script directory
rm(list=ls())

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/06_build_datasets_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('rlist')

####################################################################################
####################################################################################
## make list of datasets to process 
####################################################################################
####################################################################################

folder_names <- c("bbs/omiprep/omiprep_export_2026_07_20/", "f7/omiprep/omiprep_export_2026_07_20/", 
                    "f24/omiprep/omiprep_export_2026_07_20/", "f30/omiprep/omiprep_export_2026_07_20/", 
                      "fom1/omiprep/omiprep_export_2026_07_20/", "fom2/omiprep/omiprep_export_2026_07_20/", 
                        "g2/omiprep/omiprep_export_2026_07_20/", "tf3/omiprep/omiprep_export_2026_07_20/", "tf4/omiprep/omiprep_export_2026_07_20/")

# in this case we will use 'input' data and retain QC flag for users to exclude as appropriate 
data_files <- NULL
for (i in 1:length(folder_names)){
  data_files[i] <- paste0("../../../data/raw/", folder_names[i], "input/data.tsv")
}
data_files

samples_files <- NULL    
for (i in 1:length(folder_names)){
  samples_files[i] <- paste0("../../../data/raw/", folder_names[i], "input/samples.tsv")
}
samples_files

sample_sumstats_files <- NULL    
for (i in 1:length(folder_names)){
  sample_sumstats_files[i] <- paste0("../../../data/raw/", folder_names[i], "input/sample_summary.tsv")
}
sample_sumstats_files

feature_sumstats_files <- NULL    
for (i in 1:length(folder_names)){
  feature_sumstats_files[i] <- paste0("../../../data/raw/", folder_names[i], "qc/feature_summary.tsv")
}
feature_sumstats_files

biomarker_tags_files <- c("../../../data/raw/bbs/ALSPAC_ALBB-4_8_5-biomarker-tags.csv", "../../../data/raw/f7/ALSPAC_ALSPAC7-4_8_5-biomarker-tags.csv",
                         "../../../data/raw/f24/ALSPAC_ALTW-4_8_5-biomarker-tags.csv", "../../../data/raw/f30/f30-biomarker-tags.csv", 
                          "../../../data/raw/fom1/ALSPAC_ALSPACM1-4_8_5-biomarker-tags.csv", "../../../data/raw/fom2/ALSPAC_ALSPACM2-4_8_5-biomarker-tags.csv",
                            "../../../data/raw/g2/ALSPAC_ALG2-4_8_5-biomarker-tags_fixed.csv",
                              "../../../data/raw/tf3/ALSPAC_ALSPAC15-4_8_5-biomarker-tags.csv", "../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-biomarker-tags_new_IDs.csv")


####################################################################################
####################################################################################
## read in files needed for build
####################################################################################
####################################################################################

# feature labels 
feature_labels <- read.table(file="../../../data/raw/metadata/feature_metadata.txt",h=T,sep="\t")
dim(feature_labels)

# omiprep output
# input data 
data_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(data_list)){
  print(paste("reading in:", names(data_list[i])))
  data_list[[i]] <- read.table(file=data_files[i], h=T)
  print(dim(data_list[[i]]))
}

# sample metadata
sample_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(sample_list)){
  print(paste("reading in:", names(sample_list[i])))
  sample_list[[i]] <- read.table(file=samples_files[i], h=T,sep="\t", na.strings = c("", "NA"))
  print(dim(sample_list[[i]]))
}

# sample summary stats 
sample_sumstats_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(sample_sumstats_list)){
  print(paste("reading in:", names(sample_sumstats_list[i])))
  sample_sumstats_list[[i]] <- read.table(file=sample_sumstats_files[i], h=T,sep="\t", na.strings = c("", "NA"))
  print(dim(sample_sumstats_list[[i]]))
}

# feature summary stats
feature_sumstats_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(feature_sumstats_list)){
  print(paste("reading in:", names(feature_sumstats_list[i])))
  feature_sumstats_list[[i]] <- read.table(file=feature_sumstats_files[i], h=T,sep="\t", na.strings = c("", "NA"))
  print(dim(feature_sumstats_list[[i]]))
}

# biomarker tags (by sample) (not processed by omiprep)
biomarker_tags_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(biomarker_tags_list)){
  print(paste("reading in:", names(biomarker_tags_list[i])))
  biomarker_tags_list[[i]] <- read.table(file=biomarker_tags_files[i], h=T,sep=",", na.strings = c("", "NA"))
  print(dim(biomarker_tags_list[[i]]))
}

####################################################################################
####################################################################################
## build data file
####################################################################################
####################################################################################

for (i in 1:9) {
  data <- data_list[[i]]
  print(paste("processing data from:", names(data_list[i])))
  sample_stats <- sample_sumstats_list[[i]]
  sample <- sample_list[[i]]
  biomarker_tags <- biomarker_tags_list[[i]]
  # add missingness data from summary file - a useful sum stat
  data$proportion_missing <- sample_stats$missingness[match(data$sample_id,sample_stats$sample_id)]
  
  # add QC flags
  if (names(sample_list[i]) == "F30"){
    # F30 has sample flags in wide format - need to collapse first and change variable names to match other datasets
    sample$notes <- sample$Sample.notes
    sample$tags <- NA
    for (ind in 1:nrow(sample)){
      for (tag in 5:25) {
        if (sample[ind,tag]==1){
          sample[ind,29] = paste(sample[ind,29],names(sample)[tag],sep=";")
        }
      }
    }
    sample$tags <- sub("NA;", "", sample$tags) 
  }
  data$NHtags <- sample$tags[match(data$sample_id,sample$sample_id)]
  data$NHnotes <- sample$notes[match(data$sample_id,sample$sample_id)]
  data$omiprep_exclude <- sample$excluded[match(data$sample_id,sample$sample_id)]
  data$omiprep_reason_exclude <- sample$reason_excluded[match(data$sample_id,sample$sample_id)]
  
  # add flag for data presence (so that when all timepoints combined it is clear what is missing data versus sample not analysed)
  data$has_sample_at <- "yes"
  
  # reorder so metadata vars at the front
  data <- data |> relocate(sample_id, has_sample_at, omiprep_exclude, omiprep_reason_exclude, proportion_missing, NHtags, NHnotes)
  
  # add in biomarker tags where they are not empty for all samples
  # first restrict biomarker tag file to those with data 
  biomarker_tags <- biomarker_tags[biomarker_tags$sample_code %in% data$sample_id,]
  not_all_na <- function(x) {!all(is.na(x))}
  biomarker_tags_filtered <- biomarker_tags %>% select_if(not_all_na)
  dim(biomarker_tags_filtered)
  names(biomarker_tags_filtered)[2:ncol(biomarker_tags_filtered)] <- paste0(names(biomarker_tags_filtered)[2:ncol(biomarker_tags_filtered)],"_NHtag")
  print(dim(data))
  data <- merge(data, biomarker_tags_filtered,by.x = "sample_id", by.y = "sample_code", all.x=TRUE)
  print(dim(data))
  data_list[[i]] <- data
}

####################################################################################
####################################################################################
## prep feature metadata file for sharing with datanote & to make var labels
####################################################################################
####################################################################################

### create variable labels for features to annotate in main release dataset
# var name will be short label with clinic indicator as _CLINIC attached afterwards
# var label will be full name from feature worksheet

# populate line for glycerol in feature metadata (not in annotation file used originally - from F30 - but does appear in BBS dataset)
glycerol <- c("Glycerol", "Glycerol", "Glycerol", "mmol/l", "Glycolysis related metabolites", "NA", "Glycerol", "Glycerol (mmol/l)",
                "Glycolysis related metabolites", "NA", "Glycerol (mmol/l)", "Glycerol", "no")
feature_labels <- rbind(feature_labels,glycerol)
dim(feature_labels)


# make standard varnames and labels
var_labels_list <- list(BBS=NULL, F7 = NULL, F24 = NULL, F30 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)

for (i in 1:9) {
  print(paste("processing data from:", names(data_list[i])))
  data <- data_list[[i]]
  dataset_name <- names(data_list[i])
  feature_sumstats <- feature_sumstats_list[[i]]
  n_features <- nrow(feature_sumstats)
  var_labels <- data.frame(matrix(nrow = ncol(data), ncol = 2)) 
  names(var_labels)[1] <- "varname"
  names(var_labels)[2] <- "varlabel"
  
  # bring across varnames
  var_labels$varname <- names(data)
  # add labels for metabolite data columns
  var_labels$varlabel <- feature_labels$label[match(var_labels$varname,feature_labels$CSV.column.name)]
  # add labels for sample metadata columns
  var_labels$varlabel[1:7] <- c("sample ID (for linking to aln/qlet)", "indicator for whether a sample was analysed for the individual at the timepoint", 
                                "indicator variable for whether sample was flagged for exclusion by omiprep during qc",
                                "reason why sample was flagged for exclusion by omiprep during qc",
                                "sample proportion missingness (NA)", "sample QC tags from Nightingale Health", "sample notes from Nightingale Health")
  # add labels for biomarker tag columns
  n_tags <- nrow(var_labels)-(n_features+7)
  n_tags_label <- rep("metabolite specific QC tag from Nightingale Health",n_tags)
  var_labels$varlabel[(n_features+8):nrow(var_labels)] <- n_tags_label
  
  ### update varnames and labels in files to reflect timepoint
  names(data)[2:ncol(data)] <- paste0(names((data)[2:ncol(data)]),"_",dataset_name)
  var_labels$varname[2:nrow(var_labels)] <- paste0(var_labels$varname[2:nrow(var_labels)],"_", dataset_name)
  var_labels$varlabel[2:nrow(var_labels)] <- paste0(var_labels$varlabel[2:nrow(var_labels)],": ", dataset_name)
  feature_sumstats$varname <- paste0(feature_sumstats$feature_id, "_", dataset_name)
  
  ## check feature names map to labels and data
  print("No. of metabs in feature sumstats also in variable label file:")
  print(length(which(feature_sumstats$varname %in% var_labels$varname)))
  print("No. of metabs in feature sumstats not in variable label file:")
  print(length(which(!feature_sumstats$varname %in% var_labels$varname)))
  
  print("No. of metabs in feature sumstats also in data file:")
  print(length(which(feature_sumstats$varname %in% names(data))))
  print("No. of metabs in feature sumstats not in data file:")
  print(length(which(!feature_sumstats$varname %in% names(data))))
  
  # move dataset varname to front of sumstats
  feature_sumstats <- feature_sumstats %>%
    relocate(varname)
  
  # store
  var_labels_list[[i]] <- var_labels
  feature_sumstats_list[[i]] <- feature_sumstats
  data_list[[i]] <- data
}

####################################################################################
####################################################################################
## sample summaries 
####################################################################################
####################################################################################

# summarise data
for (i in 1:9) {
  data <- data_list[[i]]
  print(paste("processing data from:", names(data_list[i])))
  print("Number of samples:")
  print(nrow(data))
  print("Number of features:")
  n_biomarker_tags <- length(grep("_NHtag", names(data)))
  print(ncol(data)-7-n_biomarker_tags)
  print("Number of features with QC tags")
  print(n_biomarker_tags)
  print("Sample QC flags:")
  print(table(data$NHtags, useNA = "always"))
  print("Sample notes:")
  print(table(data$NHnotes))
  print("Sample QC flags & notes:")
  print(table(data$NHtags, data$NHnotes, useNA = "always"))
  print("omiprep recommended exclusions")
  print(table(data$omiprep_exclude, data$omiprep_reason_exclude, useNA = "always"))
  print("omiprep recommended exclusions & NH notes")
  print(table(data$omiprep_reason_exclude, data$NHnotes, useNA = "always"))
  print("omiprep recommended exclusions & NH tags")
  print(table(data$omiprep_reason_exclude, data$NHtags, useNA = "always"))
}

####################################################################################
####################################################################################
## save out into alspac_release dir
####################################################################################
####################################################################################

for (i in 1:9) {
  print(paste("Saving out data files for:", names(data_list[i])))
  dataset_name <- names(data_list[i])
   
  ## save data into alspac_release dir 
  write.table(data_list[[i]] , file=paste0("../../../data/alspac_release/labID/nmr_data_", dataset_name, ".txt"), sep="\t", quote=F, row.names = F)
  
  ## save feature metadata into alspac_release dir 
  write.table(feature_sumstats_list[[i]], file=paste0("../../../data/alspac_release/labID/nmr_feature_metadata_", dataset_name, ".txt"), sep="\t", quote=F, row.names = F)
  
  ## save labels into alspac_release dir 
  write.table(var_labels_list[[i]], file=paste0("../../../data/alspac_release/labID/nmr_varlabels_", dataset_name, ".txt"), sep="\t", quote=F, row.names = F)
}

# save out
list.save(data_list,"../../../data/intermediate/06_data_list.rds")
list.save(feature_sumstats_list,"../../../data/intermediate/06_feature_sumstats_list.rds")
list.save(var_labels_list,"../../../data/intermediate/06_var_labels_list.rds")

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
