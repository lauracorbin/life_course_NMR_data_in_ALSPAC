# script to fix G2 IDs prior to running metaboprep

####################
###### SET UP ######
####################

# load libraries
library('tidyr')
library('tidyverse')
library('readxl')
library('readr')
library('dplyr')

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/02_fix_g2_ids_", date, ".out"),split=TRUE)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# Read in data file in which IDs are fixed
metabdata <- read_tsv("../../../data/raw/g2/ALSPAC_ALG2-4_8_5 id fix.tsv", col_names = TRUE)

# read in metadata file in which IDs are not fixed but samples in the same order as in updated (fixed) version of data
sample_tags <- read_csv("../../../data/raw/g2/ALSPAC_ALG2-4_8_5-sample-tags.csv")

# read in biomarker (by sample) QC flags 
biomarker_tags <- read_csv("../../../data/raw/g2/ALSPAC_ALG2-4_8_5-biomarker-tags.csv")

####################################################################################
####################################################################################
## check structure of files
####################################################################################
####################################################################################

# View first few rows to verify structure
head(metabdata)
head(sample_tags)
head(biomarker_tags)

# Check structure
cat("metabdata:\n")
print(dim(metabdata))
print(names(metabdata))
print(length(unique(metabdata$sample_code)))

cat("\nSample tags file:\n")
print(dim(sample_tags))
print(names(sample_tags))
print(length(unique(biomarker_tags$sample_code)))

cat("\nBiomarker tags file:\n")
print(dim(biomarker_tags))
print(names(biomarker_tags))
print(length(unique(biomarker_tags$sample_code)))

####################################################################################
####################################################################################
## update IDs in metadata file
####################################################################################
####################################################################################

# Assume:
# - `metabdaat` contains the corrected IDs in column 1
# - `sample_tags` has one row per ID, and the ID is in the first column
# - samples in the same order in both files (as confirmed by Nightingale)

# Replace the first column of sample_tags with the fixed IDs
if (nrow(sample_tags) == nrow(metabdata)) {
  sample_tags[[1]] <- metabdata[[1]]
} else {
  stop("Number of IDs in the two files does not match.")
}

# Save the updated sample_tags
write_csv(sample_tags, "../../../data/raw/g2/ALSPAC_ALG2-4_8_5-sample-tags_fixed.csv")

# Replace the first column of biomarker_tags with the fixed IDs
if (nrow(biomarker_tags) == nrow(metabdata)) {
  biomarker_tags[[1]] <- metabdata[[1]]
} else {
  stop("Number of IDs in the two files does not match.")
}

message("Replacement complete. Output saved as 'ALSPAC_ALG2-4_8_5-sample-tags_fixed.csv'.")

# Save the updated biomarker_tags
write_csv(biomarker_tags, "../../../data/raw/g2/ALSPAC_ALG2-4_8_5-biomarker-tags_fixed.csv")

####################################################################################
####################################################################################
## check all now have the right number of unique IDs
####################################################################################
####################################################################################

cat("metabdata:\n")
print(dim(metabdata))
print(names(metabdata))
print("number of unique sample ids (should be equal to no. of rows):")
print(length(unique(metabdata$sample_code)))

cat("\nSample tags file:\n")
print(dim(sample_tags))
print(names(sample_tags))
print("number of unique sample ids (should be equal to no. of rows):")
print(length(unique(biomarker_tags$sample_code)))

cat("\nBiomarker tags file:\n")
print(dim(biomarker_tags))
print(names(biomarker_tags))
print("number of unique sample ids (should be equal to no. of rows):")
print(length(unique(biomarker_tags$sample_code)))

#################
###### END ######
#################

print(sessionInfo())

# save output
sink()

# remove data
rm(list = ls())