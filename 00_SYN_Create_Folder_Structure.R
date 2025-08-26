# This script creates the folder structure in the MODEL-AD RNA Seq Harmonization
# Project on Synapse that Nextflow and the other scripts in this project rely
# on. It only needs to be run once, unless a completely new study is added to
# the data. In that case, existing folders will remain as-is and the new study's
# folders will be added to the existing structure.

library(synapser)

synLogin(silent = TRUE)

# Set up -----------------------------------------------------------------------

# The top-level "Staging" folder in the harmonization project
data_folder_syn_id <- "syn25882599"

# The following folders were pre-created in this project by data curators:
folder_syn_ids <- list(
  metadata = "syn61850200",
  gene_expression_main = "syn64314221",
  gene_counts = "syn51132850",
  gene_norm_counts = "syn51132852"
)

# We need to add these folders to the structure:
folders_add <- list(
  bam_files = "Gene Expression (BAM files)",
  nf_pipeline = "Nextflow Pipeline Input",
  ref_genomes = "Reference Genomes"
)

# List of studies that should get folders added in various places
study_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"),
                       comment.char = "#")


# Add main folders to Data/ ----------------------------------------------------

# This will do nothing for folders that already exist except get their Synapse IDs
create_folder <- function(folder_name, parent_id) {
  folder <- Folder(name = folder_name, parent = parent_id)
  folder <- synStore(folder)
  return(folder$id)
}

folders_add_syn_ids <- list(
  bam_files = create_folder(folders_add$bam_files, folder_syn_ids$gene_expression_main),
  nf_pipeline = create_folder(folders_add$nf_pipeline, data_folder_syn_id),
  ref_genomes = create_folder(folders_add$ref_genomes, data_folder_syn_id)
)

# Combine all folders into one list
folder_syn_ids <- c(folder_syn_ids, folders_add_syn_ids)

# The Nextflow folder has non-study-specific subfolders
config <- Folder("Configuration", parent = folder_syn_ids$nf_pipeline)
sheets <- Folder("Sample Sheets", parent = folder_syn_ids$nf_pipeline)
synStore(config)
synStore(sheets)

# Currently, the metadata and reference genome folders do not have any sub-folders


# Add sub-folders for each study where needed ----------------------------------

for (study_name in study_list$Study) {
  # Folder for BAM files
  bam <- Folder(study_name, parent = folder_syn_ids$bam_files)
  bam <- synStore(bam)

  # Folder + sub-folders for raw counts files
  counts <- Folder(study_name, parent = folder_syn_ids$gene_counts)
  counts <- synStore(counts)
  qc <- Folder("Quality Control", parent = counts)
  qc <- synStore(qc)

  # TODO folders for normalized counts have not been determined yet
  # TODO folders for differential expression have not been determined yet
}
