# code to build feature metadata to go with dataset release

####################
###### SET UP ######
####################

# run code interactively from script directory

# capture output
rm(list=ls())
date <- Sys.Date()
sink(file=paste0("../../../results/11_build_feature_metadata_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('dplyr')
library('data.table')
library('openxlsx')

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in data files
data_file_names <- list.files(path = "../../../data/alspac_release/labID", pattern="metadata", full.names=TRUE)
data_file_names

data_files_list <- list(BBS=NULL, F24 = NULL, F30 = NULL, F7 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)
for (i in 1:length(data_file_names)){
  data_files_list[[i]] <- as.data.frame(fread(file=data_file_names[[i]], h=T, sep="\t"))
}

# feature labels 
feature_labels <- read.table(file="../../../data/raw/metadata/feature_metadata.txt",h=T,sep="\t")
dim(feature_labels)

####################################################################################
####################################################################################
## process
####################################################################################
####################################################################################

## populate line for glycerol in feature metadata (not in annotation file used originally - from F30 - but does appear in BBS dataset)
glycerol <- c("Glycerol", "Glycerol", "Glycerol", "mmol/l", "Glycolysis related metabolites", "NA", "Glycerol", "Glycerol (mmol/l)",
              "Glycolysis related metabolites", "NA", "Glycerol (mmol/l)", "Glycerol", "no")
feature_labels <- rbind(feature_labels,glycerol)
dim(feature_labels)

## merge metadata into summary files
metadata_list <- list(BBS=NULL, F24 = NULL, F30 = NULL, F7 = NULL, FOM1 = NULL, FOM2 = NULL, G2 = NULL, TF3 = NULL, TF4 = NULL)

for (i in 1:length(data_files_list)){
  print(paste("processing:", names(data_files_list[i])))
  # merge with feature labels
  data <- data_files_list[[i]]
  print(dim(data))
  metadata_list[[i]] <- right_join(feature_labels,data,by=c("CSV.column.name" = "feature_id")) # keep only metabolites in data file 
  print(dim(metadata_list[[i]]))
  metadata_list[[i]] <- metadata_list[[i]] %>% relocate(varname)
  # drop repeated columns
  metadata_list[[i]] <- metadata_list[[i]] %>% select(-Group, -Subgroup, -raw.label)
}

names(metadata_list[[1]])

## prep front sheet of excel file - contains labels/explainers for columns
labels_sheet <- as.data.frame(matrix(nrow=ncol(metadata_list[[1]]), ncol = 2))

labels_sheet[,1] <- names(metadata_list[[1]])
names(labels_sheet)[1] <- "column.name"
names(labels_sheet)[2] <- "column.description"

labels_sheet$column.description <- c("variable name with timepoint coding - as in data file", "metabolite name to match excel format (provided by NH)", "metabolite name to match csv format (provided by NH)",
                                     "metabolite name in full (provided by NH)", "units of metabolite measurement (provided by NH)", "variable name without timepoint coding",
                                     "class metabolite belongs to (provided by NH)","subclass metabolite belongs to (provided by NH)", "metabolite label with units included (provided by NH)",
                                     "metabolite label without units included (provided by NH)", "indicator of derived status (yes/no) (provided by NH)", 
                                     "Proportion of missing (NA) data (omiprep output)","Number of outlying values (+/- 50 IQR from median) (omiprep output)",
                                     "Sample size (non-NA) (omiprep output)", "Mean (omiprep output)", "Standard deviation (omiprep output)", 
                                     "Median (omiprep output)","Minimum value (omiprep output)", "Maximum value (omiprep output)", "Range (max - min) (omiprep output)",
                                     "Skewness (omiprep output)", "Kurtosis (omiprep output)", "Standard error (omiprep output)", "Count of NA (omiprep output)",
                                     "Variance (omiprep output)", "Variance to mean ratio (omiprep output)", "Coefficient of variation (omiprep output)",
                                     "Shaprio's W-statistic of normality on raw distribution (omiprep output)", "Shaprio's W-statistic of normality on log10 distribution (omiprep output)",
                                     "Cluster assignment (omiprep output)", "Feature is a representative independent feature (1=yes, 0=no) (omiprep output)")

# add note to say all samples are included in the summaries (so before woc and duplicate removals)
# blank row
labels_sheet <- labels_sheet %>% 
  add_row(column.name = "", column.description = "")
# note
labels_sheet <- labels_sheet %>% 
  add_row(column.name = "", 
          column.description = "PLEASE NOTE THE SUMMARY STATISTICS PRESENTED ARE EXTRACTED FROM THE DATA SUMMARY PERFORMED ON ALL ANALYSED SAMPLES (PRIOR TO QC AND OTHER EXCLUSIONS). THEREFORE, THE SAMPLE NUMBERS AND STATISTICS WILL NOT EXACTLY MATCH THE RELEASED DATA")


####################################################################################
####################################################################################
## write out 
####################################################################################
####################################################################################

# write out as two excel files - one for G0 and one for G1

g1 <- list("Column.annotations" = labels_sheet, "Before.breakfast.study" = metadata_list$BBS, "F7" = metadata_list$F7, "TF3" = metadata_list$TF3, 
                "TF4" = metadata_list$TF4, "F24" = metadata_list$F24, "F30" = metadata_list$F30)

write.xlsx(g1, file=paste0("../../../data/alspac_release/g1_feature_metadata_", date, ".xlsx"))


g0 <- list("Column.annotations" = labels_sheet, "FOM1" = metadata_list$FOM1, "FOM2" = metadata_list$FOM2)

write.xlsx(g0, file=paste0("../../../data/alspac_release/g0_feature_metadata_", date, ".xlsx"))


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
