### Find the original code here: /data/chamal/projects/yohan/projects/2022-12-transcriptome-filtering/code/workflow/get_PNC_subjects.R
### Run this code in /data/chamal/projects/yuhan/projects/2023-07-thalamus-cortex-covariance


library(tidyverse)
library(glue)

group_name <- "8_to_11"

# Get the data
behav_dict <- read_csv("resources/PNC_phenotypes/Data_Dictionary_20171108.csv")
df_qc <- read_csv("resources/PNC/PNC_ids_all_imaging_qced.csv")
df_agesex <- read_csv("resources/PNC_phenotype_data/PNC_Demographics_1598_id-sorted_phen-sex.csv") %>%
  mutate(Subject_ID=paste("PNC", SUBJID, sep="_")) %>%
  select(-SUBJID)
df_behav <- read_tsv("resources/PNC_phenotypes/PNC_phenotypes_20171031.csv") %>%
  select(dbGaP_Subject_ID:Med_Rating) %>%
  mutate(Subject_ID=paste("PNC", SUBJID, sep="_")) %>%
  select(-SUBJID)

# Join the data
df <- df_agesex %>%
  left_join(df_qc, by="Subject_ID") %>%
  left_join(df_behav, by="Subject_ID")

# Do the filtering
df_filtered <- df %>%
  filter(Med_Rating %in% c(0,1), DBMQC_Pass, CivetQC_Pass, Age_header >= 8.0, Age_header <= 11.0, Sex_phenotypes == "male")

# Output
dir.create(glue("data/PNC/populations/healthy/groups_6/male/{group_name}"))
df_filtered %>% write_csv(glue("data/PNC/populations/healthy/groups_6/male/{group_name}.csv"))
df_filtered

# Symlink subject dirs
linkdir <- glue("data/PNC/populations/healthy/groups_6/male/{group_name}/civet_outputs")
dir.create(linkdir)

for (i in 1:nrow(df_filtered)) {
  sid <- df_filtered$Subject_ID[[i]]
  target <- glue("resources/PNC/civet_output_consolidated/{sid}")
  link <- glue("{linkdir}/{sid}")
  R.utils::createLink(link=link, target=target)
}

