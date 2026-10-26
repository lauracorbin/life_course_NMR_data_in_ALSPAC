# script to deconstruct at 30 excel file into flat text format for input into metaboprep v1

# load libraries
library('tidyr')
library('tidyverse')
library('readxl')
library('readr')
library('dplyr')

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/03_extract_f30_", date, ".out"),split=TRUE)


####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in data 
results <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Results", col_names=FALSE, range = "C16:IR2909", na=c("NA", "TAG")))
results_headers <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Results", col_names=FALSE, range = "C11:IR11", col_types = "text"))

# read in tags per biomarker
tagsperbiomarker <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Tags per biomarker", col_names=FALSE, range = "B16:IQ2909", na="", col_types = "text"))
tagsperbiomarker_headers <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Tags per biomarker", col_names=FALSE, range = "B11:IQ11", na="", col_types = "text"))

# read in sample qc
sampleqc <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Quality control tags and notes", col_names=FALSE, range = "B16:Z2909", na="", col_types = "text"))
samplesqc_header <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Quality control tags and notes", col_names=FALSE, range = "B12:Z12", na="", col_types = "text"))

# read in feature annotations
feature_info <- as.data.frame(readxl::read_xlsx("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet="Biomarker annotations", col_names=TRUE, col_types = "text"))

####################################################################################
####################################################################################
## reformat to look like other timepoint datasets
####################################################################################
####################################################################################

# replace metabolite names with csv version (to match other timepoints)

# results file
featurenames_r <- as.data.frame(t(results_headers[,2:ncol(results_headers)]))
featurenames_r$csv_names <- feature_info$`CSV column name`[match(featurenames_r$V1,feature_info$`Excel column name`)]
names(results) <- c("sample_code", featurenames_r$csv_names)

# tags
featurenames_t <- as.data.frame(t(tagsperbiomarker_headers[,2:ncol(tagsperbiomarker_headers)]))
featurenames_t$csv_names <- feature_info$`CSV column name`[match(featurenames_t$V1,feature_info$`Excel column name`)]
names(tagsperbiomarker) <- c("sample_code", featurenames_t$csv_names)

# sample qc
names(sampleqc) <- c("sample_code", samplesqc_header[2:ncol(samplesqc_header)])

####################################################################################
####################################################################################
## write out
####################################################################################
####################################################################################

write.table(results, "../../../data/raw/f30/f30-data.csv",quote=F,row.names=F,sep=",")
write.table(tagsperbiomarker, "../../../data/raw/f30/f30-biomarker-tags.csv",quote=F,row.names=F,sep=",", na="")
write.table(sampleqc, "../../../data/raw/f30/f30-sample-tags.csv",quote=F,row.names=F,sep=",", na="")


#################
###### END ######
#################

print(sessionInfo())

# save output
sink()

# remove data
rm(list = ls())
