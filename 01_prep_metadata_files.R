# code to reformat feature metadata file so it resembles a post ng_anno merge datafile.

####################
###### SET UP ######
####################

# load libraries
library('tidyr')
library('tidyverse')
library('readxl')

# capture output
date <- Sys.Date()
sink(file=paste0("../../../results/01_prep_metadata_files_", date, ".out"),split=TRUE)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in metabolite metadata from F30 release (as not present for other releases)
featuredata <- as.data.frame(readxl::read_excel("../../../data/raw/f30/24701 ALSPAC30yr-21-Nov-2024-Results.xlsx", sheet = 4, na = c("")))

####################################################################################
####################################################################################
## reformat
####################################################################################
####################################################################################

# add in columns that would be added on merging with ng_anno (because we need to take that step out)
featuredata$metabolite <- featuredata$`CSV column name`
featuredata$raw.label <- paste0(featuredata$`Biomarker name`," (", featuredata$Unit, ")")
featuredata$class <- featuredata$Group
featuredata$subclass <- featuredata$Subgroup
featuredata$label <- featuredata$raw.label
featuredata$label.no.units <- featuredata$`Biomarker name`

# add derived features - use units to derive
featuredata$derived_features <- "no"
featuredata <- featuredata %>% mutate(`derived_features` = replace(`derived_features`, `Unit`=="%" | `Unit`=="ratio", "yes")) # class as single carrier

####################################################################################
####################################################################################
## save out
####################################################################################
####################################################################################

write.table(featuredata, "../../../data/raw/metadata/feature_metadata.txt", sep="\t", row.names = F, quote = F)

#################
###### END ######
#################

print(sessionInfo())

# save output
sink()

# remove data
rm(list = ls())



