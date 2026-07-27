# This script creates the folder structure in the MODEL-AD RNA Seq Harmonization
# Project on Synapse that Nextflow and the other scripts in this project rely
# on. Specifically, folders are created in Staging that mirror the folders in
# Data, and a folder for each study is added to the raw counts, bam
# files, and genotype validation folders, using the same name for the study as
# appears in annotations on Synapse.
#
# The script then compiles the IDs of all of the folders in Data, plus all of the
# mirrored folders it made, and writes them to a config file that can be used
# across the pipeline.
#
# This script only needs to be run once, unless a completely new study is added
# to the data or the overall folder structure changes. In that case, existing
# folders will remain as-is and the new study's folders will be added to the
# existing structure.
#
# The locations of the main folders are read from the config.yml file.

library(synapser)
library(stringr)
library(dplyr)

synLogin(silent = TRUE)

# Set up -----------------------------------------------------------------------

top_level_syn_ids <- config::get("top_level_syn_ids")

# List of studies that should get folders added in various places.
study_list <- config::get("studies")

# Which studies to create folders for. Can add all studies in study_list, a
# subset, only 1 study, or no studies. No studies (empty vector) just mirrors
# the data folder structure into staging without adding any study-level folders.
studies_to_add <- c()

# Helper function -- this will do nothing for folders that already exist except
# get their Synapse IDs
create_folder <- function(folder_name, parent_id) {
  folder <- Folder(name = folder_name, parent = parent_id)
  folder <- synStore(folder)
  return(folder)
}

# Mirror the folder structure for released data into Staging -------------------

## Mirror helper function ------------------------------------------------------

# This is a recursive function that moves through the sub-folders of Data/ on
# Synapse and creates folders with identical name and structure in Staging/.
# Along the way, it creates data frames with information about each folder it
# creates, and binds them all together at the end.
#
# This is significantly faster than using synapserutils::walk() or
# synapserutils::copy() because we skip any folders that are named after studies
# (and their subfolders) and don't need to filter them out / delete them on
# Synapse afterward.
#
# To make the code clearer, variable names for folders that exist in "Data" have
# `data_` in front, while those that exist in "Staging" have `staging_`.
#
# Arguments:
#   data_folder - a synapser Folder object
#   staging_parent_id - a string, the parent folder in Staging that should
#     mirror the structure of `data_folder` in Data
#   staging_folder_path - a filepath-like string used strictly for printing,
#     which gives the location of `staging_parent_id` relative to "Staging/" on
#     Synapse
#
# Returns:
#   a data frame with info about each folder we created in Staging/ and its
#   corresponding folder in Data/
mirror <- function(data_folder, staging_parent_id, staging_folder_path) {
  data_subfolders <- synGetChildren(data_folder, includeTypes = list("folder")) |>
    as.list()

  # If no folders exist inside this one, return. If any of the child folders are
  # named after studies, we've reached the lowest level we want to mirror in
  # this folder. Don't mirror the study folders, instead just return. We also
  # don't want to mirror any folders in `Custom Genome Benchmarking/Raw Counts`.
  if (length(data_subfolders) == 0 ||
      any(sapply(data_subfolders, "[[", "name") %in% study_list) ||
      grepl("Custom Genome Benchmarking\\/Raw Counts", staging_folder_path)) {
    return(NULL)
  }

  # Otherwise, continue mirroring the sub-folders
  staging_subfolders <- lapply(data_subfolders, function(data_child) {

    # New folder in Staging with the same name as `data_child`
    staging_new_folder <- create_folder(data_child$name, staging_parent_id)
    staging_new_path <- paste(staging_folder_path, staging_new_folder$name, sep = "/")

    staging_folder_df <- data.frame(name = staging_new_folder$name,
                                    id = staging_new_folder$id,
                                    parent = staging_new_folder$parentId,
                                    path = staging_new_path,
                                    # Save the original ID too
                                    released_id = data_child$id)

    # printing
    spaces <- rep(" ", 4*str_count(staging_folder_path, "/")) |> paste(collapse = "")
    print(str_glue("{spaces}|-- {staging_new_folder$id}: {staging_new_path}"))

    # Recursively mirror this sub-folder's children, if they exist
    new_children <- mirror(data_child,
                           staging_parent_id = staging_new_folder$id,
                           staging_folder_path = staging_new_path)

    return(rbind(staging_folder_df, new_children))
  })

  return(do.call(rbind, staging_subfolders))
}


## Mirror the folders ----------------------------------------------------------

data_folder <- synGet(top_level_syn_ids$released_data)
all_folders <- mirror(data_folder, top_level_syn_ids$staging, "Staging")

# Remove genome benchmarking subfolders, add a config-friendly name, add the
# path to the released data, edit some of the config names
all_folders <- subset(all_folders,
                      !grepl("Custom Genome Benchmarking\\/", path)) |>
  mutate(
    # Remove anything between parentheses, lower-case, and replace spaces with
    # underscores
    config_name = str_replace_all(name, " \\(.*\\)", "") |>
           tolower() |>
           str_replace_all(" ", "_"),
    # Shorten some of the names and add nf_ in front of Nextflow fields for clarity
    config_name = case_match(
      config_name,
      "differential_expression_analysis" ~ "de_analysis",
      "nextflow_pipeline_input" ~ "nf_pipeline_input",
      "configuration" ~ "nf_configuration",
      "sample_sheets" ~ "nf_samplesheets",
      .default = config_name
    ),
    # Path in released Data/ folder is a straight replacement of Staging => Data
    # in all the file paths
    released_path = str_replace(path, "Staging", "Data")
  )

# Write to a config file. We want to include comments, so we write manually
# instead of using write_yaml. Not all the folders listed are used in the
# pipeline but it's simpler to write everything instead of manually creating a
# list. Each line is padded with whitespace so syanpse IDs appear nicely aligned
# in the file.
# File structure is as follows:
# ```
# # Auto-generated file
# default:
#
#   # IDs for staging folders (upload or download)
#   staging_syn_ids:
#     <name>: "<id>"  # <path>
#     ...
#
#   # IDs for released data folders (download only)
#   released_data_syn_ids:
#     <name>: "<id>"  # <path>
#     ...
# ```
pad <- function(config_name, id, path) {
  base_string <- str_glue("    {config_name}:") |>
    str_pad(24, side = "right")
  str_glue("{base_string} \"{id}\"  # {path}")
}

staging <- sapply(1:nrow(all_folders), function(N) {
  row <- all_folders[N, ]
  pad(row$config_name, row$id, row$path)
})

# Manually add a norm_counts staging folder, which is never released to Data/
# and is only used to temporarily store norm counts for the explorer
norm_counts <- create_folder("Explorer Normalized Counts",
                             parent_id = top_level_syn_ids$staging)

staging <- c(
  staging,
  pad("norm_gene_counts", norm_counts$id, str_glue("Staging/{norm_counts$name}"))
)

# Data folder IDs
data <- sapply(1:nrow(all_folders), function(N) {
  row <- all_folders[N, ]
  pad(row$config_name, row$released_id, row$released_path)
})

lines <- c("# Auto-generated file",
           "default:", "",
           "  # IDs for staging folders (upload or download)",
           "  staging_syn_ids:", staging, "",
           "  # IDs for released data folders (download only)",
           "  released_data_syn_ids:", data)

cfg_file <- file("config_project_syn_ids.yml", "w")
writeLines(lines, con = cfg_file)
close(cfg_file)


# Add sub-folders for each study where needed ----------------------------------

# The new IDs should be auto-pulled into the main config now
staging_ids <- config::get("staging_syn_ids")

for (study_name in studies_to_add) {
  # Folder for BAM files
  bams <- create_folder(study_name, parent_id = staging_ids$bam_files)

  # Folder + sub-folders for raw counts files
  counts <- create_folder(study_name, parent_id = staging_ids$raw_gene_counts)
  qc <- create_folder(study_name, parent_id = staging_ids$quality_control)

  # Folder for genotype validation
  geno <- create_folder(study_name, parent_id = staging_ids$genotype_validation)
}
