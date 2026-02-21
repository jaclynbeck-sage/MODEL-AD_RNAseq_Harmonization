# This script creates the folder structure in the MODEL-AD RNA Seq Harmonization
# Project on Synapse that Nextflow and the other scripts in this project rely
# on. Specifically, a folder for each study is added to the raw counts, bam
# files, and genotype validation folders, using the same name for the study as
# appears in annotations on Synapse.
#
# This script only needs to be run once, unless a completely new study is added
# to the data. In that case, existing folders will remain as-is and the new
# study's folders will be added to the existing structure.
#
# The locations of the main folders are read from the config.yml file.

library(synapser)

synLogin(silent = TRUE)

# Set up -----------------------------------------------------------------------

folder_syn_ids <- config::get("folder_syn_ids", config = "default")

# List of studies that should get folders added in various places
study_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"),
                       comment.char = "#")

# Helper function -- this will do nothing for folders that already exist except
# get their Synapse IDs
create_folder <- function(folder_name, parent_id) {
  folder <- Folder(name = folder_name, parent = parent_id)
  folder <- synStore(folder)
  return(folder$id)
}


# Add sub-folders for each study where needed ----------------------------------

for (study_name in study_list$Study) {
  # Folder for BAM files
  bams <- create_folder(study_name, parent_id = folder_syn_ids$bams)

  # Folder + sub-folders for raw counts files
  counts <- create_folder(study_name, parent_id = folder_syn_ids$raw_counts)
  qc <- create_folder("Quality Control", parent_id = counts)

  # Folder for genotype validation
  geno <- create_folder(study_name, parent_id = folder_syn_ids$genotype_validation)
}
