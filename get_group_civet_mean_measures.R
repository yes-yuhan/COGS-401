#####################
# This script gets the mean CIVET measures for a specific population group
#####################

library(glue)

sex <- "male"
age_group <- "10_to_13"

path_civet <- "resources/civet_mean_measures.csv"

path_group <- glue(paste0("analyses/PNC/populations/healthy/groups/{sex}/{age_group}", ".csv"))

civet_measure_means_file <- read.csv(path_civet)
population_group_file <- read.csv(path_group)

group_IDs <- subset(population_group_file, select = "Subject_ID")
ages <- subset(population_group_file, select = "Age_header")

group_civet <- merge(civet_measure_means_file, group_IDs, by = "Subject_ID")

group_civet_age <- cbind(group_civet, ages)

write.csv(group_civet_age, glue(paste0("analyses/matrices/civet_mean_measures/{sex}/{age_group}",".csv")))


