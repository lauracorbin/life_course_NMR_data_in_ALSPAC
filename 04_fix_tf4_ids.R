# script to deal with duplicate sample IDs TF4
# one per pair removed - one with least missingness

####################
###### SET UP ######
####################

# Load necessary library
library(dplyr)

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/04_fix_tf4_ids_", date, ".out"),split=TRUE)

##########################
###### read in data ######
##########################

# Read the results file
df <- read.csv("../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5.csv", stringsAsFactors = FALSE)
head(df)
dim(df)

# read the sample tags
sample_tags <- read.csv("../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-sample-tags.csv", stringsAsFactors = FALSE)
head(sample_tags)
dim(sample_tags)

# read the biomarker tags
biom_tags <- read.csv("../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-biomarker-tags.csv", stringsAsFactors = FALSE)
head(biom_tags)
dim(biom_tags)

###############################
###### process data file ######
###############################

print("Number of unique sample IDs in data")
length(unique(df$sample_code))

# Identify duplicated sample_codes (both first and second occurrences)
duplicated_codes <- df$sample_code[duplicated(df$sample_code) | duplicated(df$sample_code, fromLast = TRUE)]

# Extract all rows with duplicated sample_codes
df_duplicates <- df %>% filter(sample_code %in% duplicated_codes)
dim(df_duplicates)

# Save the duplicate sample IDs out
write.csv(df_duplicates[,c("sample_code")], "../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5_duplicates.csv", row.names = FALSE)

# add row number to update sample IDs so no longer matching
df$ID <- seq.int(nrow(df))
df$sample_code <- paste0(df$sample_code,"_",df$ID)

# remove index column
df <- df[,c(1:(ncol(df)-1))]

# make rownames sample IDs
rownames(df) <- df$sample_code

# Save modified dataset
write.csv(df[,c(2:ncol(df))], "../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5_new_IDs.csv", row.names = TRUE)

######################################
###### process sample tags file ######
######################################

# add row number to update sample IDs so no longer matching
sample_tags$ID <- seq.int(nrow(sample_tags))
sample_tags$sample_code <- paste0(sample_tags$sample_code,"_",sample_tags$ID)

# Save cleaned
write.csv(sample_tags[,c(1:(ncol(sample_tags)-1))], "../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-sample-tags_new_IDs.csv", row.names = FALSE, quote=F)

#########################################
###### process biomarker tags file ######
#########################################

# add row number to update sample IDs so no longer matching
biom_tags$ID <- seq.int(nrow(biom_tags))
biom_tags$sample_code <- paste0(biom_tags$sample_code,"_",biom_tags$ID)

# Save cleaned
write.csv(biom_tags[,c(1:(ncol(biom_tags)-1))], "../../../data/raw/tf4/ALSPAC_ALSPAC17-4_8_5-biomarker-tags_new_IDs.csv", row.names = FALSE, quote=F)

#################
###### END ######
#################

print(sessionInfo())

# save output
sink()

# remove data
rm(list = ls())
