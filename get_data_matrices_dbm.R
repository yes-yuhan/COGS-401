#!/usr/bin/env Rscript
#
# This script grabs all the Jacobian determinants under cortical and thalamic
# masks, collects them into [subject x voxel] matrices, normalizes them by
# regressing out the total thalamic or cortical volume and age, and saves them as npy
# objects on disk so that they can be further worked on in either R (via the
# RcppCNPy library) or Python (via numpy)

library(tidyverse)
library(glue)
library(RMINC)
library(normtools)
library(RcppCNPy)
library(pbapply)

# Define age and groups
sex <- "male"
age_group <- "10_to_13"

# Inputs
atlas_tag <- "AHRA_annotation_full"
thalamus_jd_tag <- "jd_full_0mm"
cortex_jd_tag <-  "jd_full_0mm_lowres"

# Paths
thalamus_mask_file <- glue("resources/yohan/derivatives/PNC/masks/{atlas_tag}/thalamus_mask_highres.mnc")
cortex_mask_file <- glue("resources/yohan/derivatives/PNC/masks/{atlas_tag}/cortex_mask_lowres.mnc")

thalamus_jd_dir <- glue("resources/yohan/derivatives/PNC/dbm/{thalamus_jd_tag}")
cortex_jd_dir <- glue("resources/yohan/derivatives/PNC/dbm/{cortex_jd_tag}")

outfile_dat_thalamus <- glue("analyses/matrices/{atlas_tag}/{sex}/{age_group}/dat_thalamus_{thalamus_jd_tag}_normalized_tvolreg_agereg.npy")
outfile_dat_cortex <- glue("analyses/matrices/{atlas_tag}/{sex}/{age_group}/dat_cortex_{cortex_jd_tag}_normalized_cvolreg_agereg.npy")


thalamus_cortex_summaries_file <- glue(paste0("analyses/matrices/{atlas_tag}/structural_summaries/dbm_thalamus_cortex_volumes/{sex}/{age_group}",".csv"))
subjects_dbm_file <- glue(paste0("analyses/PNC/populations/healthy/groups/{sex}/{age_group}",".csv"))

####################
# Read data
####################

# Read input subjects
subjects_dbm <- read_csv(subjects_dbm_file)

# Read volumes for normalization
# and confirm subjects are in the same order
vf <- read_csv(thalamus_cortex_summaries_file)
if (!all(vf$Subject_ID==subjects_dbm$Subject_ID)) {
  stop("Error - subject ordering issue")
}

####################
# Setup output objects
####################

# Read masks
tmvol <- mincGetVolume(thalamus_mask_file)
cmvol <- mincGetVolume(cortex_mask_file)

# Get matrix sizes
n <- nrow(subjects_dbm)
p_thal <- length(which(tmvol > 0.5))
p_cortex <- length(which(cmvol > 0.5))

# Empty matrices
dat_thalamus <- matrix(nrow=n, ncol=p_thal)
dat_cortex <- matrix(nrow=n, ncol=p_cortex)
rownames(dat_thalamus) <- rownames(dat_cortex) <- subjects_dbm$Subject_ID

####################
# Read in DBM data into matrix
####################

# Read in data
# 10-15 minutes
start_time <- Sys.time()
prog <- txtProgressBar(max=n, style=3)
for (i in 1:n) {
  # This ID
  subject_id <- subjects_dbm$Subject_ID[[i]]

  # File paths
  #
  # If file system seems unstable, use the loop below instead:
  jd_thal_file <- jd_cortex_file <- character(0)
  while (!((length(jd_thal_file)==1) & (length(jd_cortex_file)==1))) {
   jd_thal_file <- list.files(thalamus_jd_dir, pattern=subject_id, full.names = T)
   jd_cortex_file <- list.files(cortex_jd_dir, pattern=subject_id, full.names = T)
  }
  #jd_thal_file <- list.files(thalamus_jd_dir, pattern=subject_id, full.names = T)
  #jd_cortex_file <- list.files(cortex_jd_dir, pattern=subject_id, full.names = T)

  # Get data
  tvec <- mincGetVolume(jd_thal_file)[tmvol > 0.5]
  cvec <- mincGetVolume(jd_cortex_file)[cmvol > 0.5]

  # Fill in matrix, line by line
  dat_thalamus[i,] <- tvec
  dat_cortex[i,] <- cvec

  # Output progress
  setTxtProgressBar(prog, i)
}
close(prog)
end_time <- Sys.time()
tdiff <- end_time - start_time
print(glue("Compute time: {tdiff}"))

# Double check ordering again
all(rownames(dat_thalamus)==vf$Subject_ID)
all(rownames(dat_cortex)==vf$Subject_ID)

####################
# Normalize
####################

# Normalize
# takes seconds to minutes
print(Sys.time())
norm_dat_thalamus <- norm_by_model(dat_thalamus, ~ Thalamus + Age_header, df = vf)
print(Sys.time())
norm_dat_cortex <- norm_by_model(dat_cortex, ~ Cortex + Age_header, df = vf)
print(Sys.time())

# Write out data
npySave(filename = outfile_dat_thalamus, object = norm_dat_thalamus)
npySave(filename = outfile_dat_cortex, object = norm_dat_cortex)

# Done
Sys.time()

