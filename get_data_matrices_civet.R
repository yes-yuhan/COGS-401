#!/usr/bin/env Rscript

library(tidyverse)
library(glue)
library(normtools)
library(RcppCNPy)
library(pbapply)
# Inputs
atlas_tag <- "AHRA_annotation_full"
measure_tags <- list(
  gyrification_left=c("surfaces", "gyrification_rsl_left", 40962),

  area_left=c("surfaces", "mid_surface_rsl_left_native_area_0mm", 40962),

  volume_left=c("surfaces", "surface_rsl_left_native_volume_0mm", 40962),

  mc_mid_left=c("thickness", "native_mc_rsl_0mm_mid_left", 40962),

  thickness_left=c("thickness", "native_rms_rsl_tlaplace_0mm_left", 40962)
 
) # Format is a vector: civet subfolder, tag, number of vertices
#  area_right=c("surfaces", "mid_surface_rsl_right_native_area_0mm", 40962),
#  asym_pos=c("surfaces", "native_pos_rsl_asym_hemi", 40962),
#  volume_right=c("surfaces", "surface_rsl_right_native_volume_0mm", 40962),
#  mc_gray_left=c("thickness", "native_mc_rsl_0mm_gray_left", 40962),
#  mc_gray_right=c("thickness", "native_mc_rsl_0mm_gray_right", 40962),
#  mc_mid_right=c("thickness", "native_mc_rsl_0mm_mid_right", 40962),
#  mc_white_left=c("thickness", "native_mc_rsl_0mm_white_left", 40962),
#  mc_white_right=c("thickness", "native_mc_rsl_0mm_white_right", 40962),
#  asym_thickness=c("thickness", "native_rms_rsl_tlaplace_0mm_asym_hemi", 40962),
#  gyrification_right=c("surfaces", "gyrification_rsl_right", 40962), 
# thickness_right=c("thickness", "native_rms_rsl_tlaplace_0mm_right", 40962)

# Define the sex and age group
sex <- "male"
age_group <- "10_to_13"


# A file that has the civet mean measures for the target subjects
civet_measure_means_file <- glue(paste0("analyses/matrices/civet_mean_measures/{sex}/{age_group}", ".csv"))

# A file that has a list of subject IDs for the included subjects
# Not sure whether to include this or not since all subjects are already in the civet mean measures file
# subjects_civet_file <- glue("path")


outdir <- glue("analyses/matrices/{atlas_tag}/{sex}/{age_group}")

####################
# Read data
####################

# Read input subjects
subjects_civet <- read_csv(civet_measure_means_file)

# Read volumes for normalization
# and confirm subjects are in the same order
vf <- read_csv(civet_measure_means_file)
#if (!all(vf$Subject_ID==subjects_civet$Subject_ID)) {
#  stop("Error - subject ordering issue")
#}
names(vf)
####################
# Loop over input measures
####################

for (m in 1:length(measure_tags)) {

  # Get info on measure
  measure <- names(measure_tags)[[m]]
  civet_output_subdir <- measure_tags[[measure]][[1]]
  measure_tag <- measure_tags[[measure]][[2]]
  n_vertices <- as.integer(measure_tags[[measure]][[3]])

  # Template for file location
  infile_template <- glue("resources/yohan/derivatives/PNC/civet_output_consolidated/{{subject_id}}/{civet_output_subdir}/{{subject_id}}_{measure_tag}.txt")

  # Output file
  outfile_measure <- glue("{outdir}/dat_cortex_{measure}_cmeanreg_agereg.npy")

  # Status
  print(glue("[{Sys.time()}] Working on measure {m} of {length(measure_tags)}: {measure}"))
  print(glue("> CIVET data location: {infile_template}"))
  print(glue("> Number of vertices: {n_vertices}\n\n"))

  # Empty matrix
  dat_measure <- matrix(nrow=nrow(subjects_civet), ncol=n_vertices)
  rownames(dat_measure) <- subjects_civet$Subject_ID

  # Read in data
  # This internal loop should take seconds to minutes
  prog <- txtProgressBar(max=nrow(subjects_civet), style=3)
  for (i in 1:nrow(subjects_civet)) {

    # This ID
    subject_id <- subjects_civet$Subject_ID[[i]]

    # File
    measure_file <- glue(infile_template)

    # Read data
    fvec <- read_csv(measure_file, col_names = "feature", col_types = "d")

    # Populate matrix
    dat_measure[i,] <- fvec$feature

    # Progress
    setTxtProgressBar(prog, i)
  }
  close(prog)

  # Check ordering again
#  if (!all(rownames(dat_measure)==vf$Subject_ID)) {
#    stop("Error - subject ordering issue")
#  }

  # Normalize
  # Also takes a few minutes
  print(glue("[{Sys.time()}] Normalizing matrix..."))
  norm_dat_measure <- norm_by_model(dat_measure, ~ Feature + Age_header, df=(vf %>% select(Subject_ID, Feature=!!sym(measure), Age_header)))

  # Write out data
  # Should take seconds
  print(glue("[{Sys.time()}] Writing out normalized data matrix...\n\n"))
  npySave(filename = outfile_measure, object = norm_dat_measure)
}

# Done
Sys.time()

