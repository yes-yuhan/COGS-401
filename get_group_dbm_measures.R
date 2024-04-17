#####################
# This script gets the mean CIVET measures for a specific population group
#####################

library(glue)

sex <- "male"
group <- "10_to_13"

path_dbm <- "resources/dbm_thalamus_cortex_volumes.csv"
path_group <- glue(paste0("analyses/PNC/populations/healthy/groups/{sex}/{group}", ".csv"))

dbm_file <- read.csv(path_dbm)
population_group_file <- read.csv(path_group)


group_IDs <- subset(population_group_file, select = "Subject_ID")
ages <- subset(population_group_file, select = "Age_header")

group_dbm <- merge(dbm_file, group_IDs, by = "Subject_ID")
group_dbm_age <- cbind(group_dbm, ages)

write.csv(group_dbm_age, glue(paste0("analyses/matrices/AHRA_annotation_full/structural_summaries/dbm_thalamus_cortex_volumes/{sex}/{group}", ".csv")))


