# code to ID duplicates and check correlations

## input = metabolite data on aln/qlet IDs
## output = correlation stats 

####################
###### SET UP ######
####################

# run code interactively from script directory

# capture output
rm(list=ls())
date <- Sys.Date()
sink(file=paste0("../../../results/08_check_dupes_", date, ".out"),split=TRUE)

# load libraries
library('tidyverse')
library('dplyr')
library('ggplot2')
library('blandr')
library('rlist')
library('irr')
library('data.table')
library('purrr')

# set seed so random selection consistent
set.seed(20260512)

####################################################################################
####################################################################################
## read in data
####################################################################################
####################################################################################

# read in dataset list from 07
data_list_orig <- readRDS("../../../data/intermediate/07_linked_dataset_list.rds")
summary(data_list_orig)

# remove G2 as this dataset has a lot of duplicates and needs to be looked at separately
# remove BBS baseline and postglu as no duplicates
data_list <- data_list_orig[c(1:5,7:8)]
summary(data_list)

# bring in feature summaries - needed to id independent features
feature_files <- readRDS("../../../data/intermediate/06_feature_sumstats_list.rds")
summary(feature_files)

# make list match order of data
feature_files <- feature_files[c(3,4,2,5,6,8,9)]
summary(feature_files)

####################################################################################
####################################################################################
## check correlations between duplicates
####################################################################################
####################################################################################

# set up dataset list for output 
corstats <- list(F24 = NULL, F30 = NULL, F7 = NULL, 
                 FOM1 = NULL, FOM2 = NULL, TF3 = NULL, TF4 = NULL)

# loop over datasets - with z-scoring
pdf("../../../results/figures/08_dups_histograms.pdf")
for (i in 1:7){
  print("processing:")
  print(names(data_list[i]))
  corstats_file <- data.frame(alnqlet = character(),sample1=character(),sample2=character(),Corstat.zscored=numeric())
  data <- as.data.frame(data_list[[i]])
  dim(data)
  feature_data <- feature_files[[i]]
  # drop those with all missing data
  target_col <- grep("omiprep_exclude_", colnames(data), value = TRUE)
  data <- data[data[[target_col]] == FALSE,]
  dim(data)
  # restrict to independent features
  tokeep <- feature_data[feature_data$independent_features==1,"feature_id"]
  tokeep <- paste0(tokeep,"_",names(data_list[i]))
  cols_tokeep <- which(colnames(data) %in% tokeep)
  data <- data[,c(1:2,cols_tokeep)]
  # scale
  data_z <- apply(data[, c(3:ncol(data))], 2, scale)
  data_z <- cbind(data[,c(1,2)],data_z)
  # id those with duplicates and check corr
  dups <- as.vector(data_z[data_z$count_iids == 2 | data_z$count_iids == 3, c("alnqlet")])
  dups_data <- data_z[data_z$alnqlet %in% dups,]
  # extract random pairs to compare
  nondups_data <- data_z[!data_z$alnqlet %in% dups,]
  # randomly select 200 pairs
  nondups_data_random <- nondups_data %>% slice_sample(n = 200)
  # for correlation - join dups and nondups
  for_cor_test <- rbind(dups_data,nondups_data_random)
  # calculate correlations
  row_correlation_matrix <- cor(t(for_cor_test[,c(3:ncol(for_cor_test))]), use="pairwise.complete.obs")
  # extract lower triangle indices
  lt_idx <- lower.tri(row_correlation_matrix, diag = FALSE)
  # Build a clean, long data frame
  df_long <- data.frame(
    Var1 = rownames(row_correlation_matrix)[row(row_correlation_matrix)[lt_idx]],
    Var2 = colnames(row_correlation_matrix)[col(row_correlation_matrix)[lt_idx]],
    Correlation = row_correlation_matrix[lt_idx]
  )
  for_cor_test$row_number <- row.names(for_cor_test)
  df_long$id1 <- for_cor_test$alnqlet[match(df_long$Var1, for_cor_test$row_number)]
  df_long$id2 <- for_cor_test$alnqlet[match(df_long$Var2, for_cor_test$row_number)]
  # identify dups
  df_long$dups <- "no"
  df_long <- df_long %>% mutate(`dups` = replace(`dups`, `id1` == `id2`, "yes"))
  # Summarize by factor
  df_long_notdups <- df_long[df_long$dups == "no",]
  print("summary of random pairwise correlation")
  print(summary(df_long_notdups$Correlation))
  hist(df_long_notdups$Correlation, main = paste0(names(data_list[i]), " - random pairs"))
  df_long_dups <- df_long[df_long$dups == "yes",]
  print("summary of duplicates pairwise correlation")
  print(summary(df_long_dups$Correlation))
  hist(df_long_dups$Correlation, main = paste0(names(data_list[i]), " - duplicate pairs"))
  corstats[[i]] <- df_long_dups[,c(3:5)] %>% 
    arrange(., Correlation) %>% 
    relocate(id1, id2, Correlation)
  corstats[[i]]$timepoint <- paste0(names(data_list[i]))
  corstats[[i]]$check_linkage <- "no"
  corstats[[i]] <- corstats[[i]] %>% mutate(`check_linkage` = replace(`check_linkage`, Correlation < 0.5, "yes"))
}  
dev.off()

####################################################################################
####################################################################################
## save out 
####################################################################################
####################################################################################

# save out cor stats
list.save(corstats,"../../../data/intermediate/08_corstats_list.rds")
merged_corstats <- bind_rows(corstats)
write.table(merged_corstats, "../../../data/intermediate/08_corstats_list_to_check.csv", sep=",", row.names = F, quote=F)

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
