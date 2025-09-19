# This script finds all RSEM output files from all studies and merges them
# together so there is one file containing all studies, per type of output data
# (gene counts, gene tpm, transcript counts, or transcript tpm).
#
# While specimen IDs within an individual study are unique, across all studies
# they are not guaranteed to be unique and there are a few specimenIDs that are
# used in two different studies. Specimen IDs in these merged files are
# therefore altered to ensure uniqueness: The study name is prepended before the
# specimen ID. The altered specimenIDs should match the "merged_file_specimenID"
# column in the merged metadata file.
#
# This script is a little fragile because of string matching and assumes:
#   1. The RSEM counts files all exist in a folder under syn51132850, and that
#      folder is named with the study exactly as specified in the merged
#      metadata file.
#   2. The RSEM counts files are all named "<study>.gene_counts.tsv",
#      "<study>.gene_tpm.tsv", "<study>.transcript_counts.tsv", or
#      "<study>.transcript_tpm.tsv"
#   3. The first two columns of each file are either ["gene_id" and
#      "transcript_id(s)"] for gene_* files, or ["transcript_id" and "gene_id"]
#      for transcript_* files.

library(synapser)
library(stringr)
library(dplyr)
library(purrr)

folder_syn_ids <- config::get("folder_syn_ids", config = "default")

syn_id_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"),
                        comment.char = "#")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url", config = "default"),
                 "/blob/main/04_SYN_MergeCountsFiles.R")
tmp_dir <- file.path("output", "tmp")

# Get all the folders of counts files that exist on Synapse
study_folders <- synGetChildren(folder_syn_ids$raw_counts,
                                includeTypes = list("folder"))$asList()
study_names <- sapply(study_folders, "[[", "name")

# Subset to specific studies here, if necessary. Currently we ignore the
# benchmarking folder.
study_folders <- study_folders[!grepl("Benchmarking", study_names)]
study_names <- sapply(study_folders, "[[", "name")

# Print a warning if there are folders on Synapse that don't correspond to the
# studies in the syn ID list file.
if (any(!(study_names %in% syn_id_list$Study))) {
  missing <- setdiff(study_names, syn_id_list$Study)
  msg <- paste0("Synapse folder(s) found for studies that do not exist in ",
                "Model_AD_SynID_list.csv: \n",
                paste(missing, collapse = ", "),
                "\nThese studies will be ignored.")
  message(msg)

  study_folders <- study_folders[study_names %in% syn_id_list$Study]
  study_names <- sapply(study_folders, "[[", "name")
}


# Loop through each study folder -----------------------------------------------

counts_list <- lapply(study_folders, function(study) {
  study_files <- synGetChildren(study$id,
                                includeTypes = list("file"))
  study_files <- study_files$asList()

  if (length(study_files) == 0) {
    message(str_glue("WARNING: {study$name} raw counts folder is empty. Skipping..."))
    return(NULL)
  }

  # Get whether each file is gene_counts, gene_tpm, transcript_counts, or
  # transcript_tpm
  file_names <- sapply(study_files, "[[", "name")
  names(study_files) <- str_extract(file_names, "(gene|transcript)_(counts|tpm)")


  ## Loop through each RSEM file -----------------------------------------------

  file_contents <- lapply(study_files, function(sf) {
    syn_file <- synGet(sf$id, downloadLocation = tmp_dir, ifcollision = "overwrite.local")
    print(syn_file$name)

    counts <- read.table(syn_file$path, header = TRUE, sep = "\t")

    # Specimen IDs are not guaranteed to be unique across all studies so we add
    # the study name to the specimen ID in the column names to ensure uniqueness
    # when all studies are merged, even though this makes the names much longer.
    # The first two columns are always gene_id and transcript_id, so we ignore
    # them
    spec_ids <- str_replace(colnames(counts)[c(-1, -2)], "^X", "")
    spec_ids <- paste(study$name, spec_ids, sep = ".")
    colnames(counts)[c(-1, -2)] <- make.names(spec_ids)

    # Temporarily change the "transcript_id.s." column in gene_counts and
    # gene_tpm files to "transcript_id" so all files have the same two
    # column names, for use below.
    if ("transcript_id.s." %in% colnames(counts)) { # gene counts
      counts <- dplyr::rename(counts, transcript_id = transcript_id.s.)
    }

    return(list("counts" = counts, "provenance" = sf$id))
  })

  # file_contents is a named list where each item corresponds to a single file,
  # and the item is a list containing "counts" and "provenance". The name of
  # each item will be either "gene_counts", "gene_tpm", "transcript_counts",
  # or "transcript_tpm".
  return(file_contents)
})

counts_list <- counts_list[lengths(counts_list) > 0] # Remove any NULL entries

# counts_list is a list where each item corresponds to a single study, and the
# item is a file_contents list as described above.


# Merge all files for each file type -------------------------------------------

# file_type is either "gene_counts", "gene_tpm", "transcript_counts", or
# "transcript_tpm". We pull the names from the first item of counts_list, but
# all items should have the same 4 names.
for (file_type in names(counts_list[[1]])) {
  # Extract all items named with <file_type> from each sub-list in counts_list
  file_info <- lapply(counts_list, "[[", file_type)
  print(file_type)

  # Extract the "counts" data frame from each sub-list in file_info, then
  # merge all the data frames together using gene_id and transcript_id as keys.
  counts <- lapply(file_info, "[[", "counts")
  counts <- purrr::reduce(counts, merge, by = c("gene_id", "transcript_id"),
                          all = TRUE)

  # Set NA values to 0: Some data sets were run before adding CLU to the genome
  # and won't have CLU in the set of genes.
  counts[is.na(counts)] <- 0

  # Extract the "provenance" vector from each sub-list in file_info
  provenance <- sapply(file_info, "[[", "provenance")

  # The first two columns (gene_id and transcript_id) may have switched position
  # during the merge operation, so this puts them back in the expected order
  # based on whether this is a gene_* file or a transcript_* file. It also
  # renames the "transcript_id" column in gene_* files back to its original name
  if (str_detect(file_type, "gene")) {
    counts <- select(counts, gene_id, transcript_id, 3:ncol(counts))
    colnames(counts)[2] <- "transcript_id(s)"
  } else if (str_detect(file_type, "transcript")) {
    counts <- select(counts, transcript_id, gene_id, 3:ncol(counts))
  }


  ## Write to file and upload to Synapse ---------------------------------------

  file_name <- file.path(tmp_dir, str_glue("Model-AD_all_studies.{file_type}.tsv"))
  write.table(counts, file_name, quote = FALSE, row.names = FALSE, sep = "\t")

  syn_file <- File(file_name, parent = folder_syn_ids$raw_counts)
  synStore(syn_file, forceVersion = FALSE, used = provenance,
           executed = github)
}


