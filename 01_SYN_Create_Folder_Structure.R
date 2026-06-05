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

project_syn_ids <- config::get("project_syn_ids")

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

# This is significantly faster than using synapserutils::walk() because we skip
# any folders that are named after studies and their subfolders, and don't need
# to parse through the structure that is returned by walk().
mirror <- function(item, staging_parent_id, folder_path) {
  children <- synGetChildren(item, includeTypes = list("folder")) |>
    as.list()

  # If no folders exist inside this one, return. If any of the child folders are
  # named after studies, we've reached the lowest level we want to mirror in
  # this folder. Don't mirror the study folders, instead just return. We also
  # don't want to mirror any folders in `Custom Genome Benchmarking/Raw Counts`.
  if (length(children) == 0 ||
      any(sapply(children, "[[", "name") %in% study_list) ||
      grepl("Custom Genome Benchmarking\\/Raw Counts", folder_path)) {
    return(NULL)
  }

  # Otherwise, continue mirroring
  sub_folders <- lapply(children, function(child) {
    new_folder <- create_folder(child$name, staging_parent_id)
    new_path <- paste(folder_path, new_folder$name, sep = "/")

    # printing
    spaces <- rep(" ", 4*str_count(folder_path, "/")) |> paste(collapse = "")
    print(str_glue("{spaces}|-- {new_folder$id}: {new_path}"))

    # Mirror children
    new_children <- mirror(child,
                           staging_parent_id = new_folder$id,
                           folder_path = new_path)

    new_folder_df <- data.frame(name = new_folder$name,
                                id = new_folder$id,
                                parent = new_folder$parentId,
                                path = new_path,
                                # Save the original ID too
                                released_id = child$id)

    return(rbind(new_folder_df, new_children))
  })

  return(do.call(rbind, sub_folders))
}

## Mirror the folders ----------------------------------------------------------

top_level_folder <- synGet(project_syn_ids$released_data)
all_folders <- mirror(top_level_folder, project_syn_ids$staging, "Staging")

# Remove genome benchmarking subfolders, add a config-friendly name, add the
# path to the released data, edit some of the config names
all_folders <- subset(all_folders, !grepl("Custom Genome Benchmarking\\/", path)) |>
  mutate(config_name = str_replace_all(name, "\\(|\\)", "") |>
           tolower() |>
           str_replace_all(" ", "_"),
         config_name = case_match(
           config_name,
           "differential_expression_analysis" ~ "de_analysis",
           "nextflow_pipeline_input" ~ "nf_pipeline_input",
           "configuration" ~ "nf_configuration",
           "sample_sheets" ~ "nf_sample_sheets",
           .default = config_name
         ),
         released_path = str_replace(path, "Staging", "Data"))

# Write to a config file. We want to include comments so we write manually
# instead of using write_yaml. Not all the folders listed are used in the
# pipeline but it's simpler to write everything instead of manually creating a
# list.
staging <- sapply(1:nrow(all_folders), function(N) {
  row <- all_folders[N, ]
  str_glue("    {row$config_name}: \"{row$id}\" \t# {row$path}")
})

data <- sapply(1:nrow(all_folders), function(N) {
  row <- all_folders[N, ]
  str_glue("    {row$config_name}: \"{row$released_id}\" \t# {row$released_path}")
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

# The new IDs should get auto-pulled into the main config now
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
