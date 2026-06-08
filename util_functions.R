library(purrr)

# Gets all child folders/files from a given directory on Synapse, fetching from
# the mirrored directories in both Staging and Data. In the case where a file
# exists in both Staging and Data, the child in Data is dropped in favor of the
# latest, unreleased version (in Staging). A message prints out when this
# happens.
#
# `config_item_name` is the name of the folder listed in
# `config_project_syn_ids.yml`, which should be the same in both the
# staging_syn_ids and released_data_syn_ids lists.
#
# This function returns a list as output by synGetChildren, with any duplicate
# items from Data removed.
syn_get_unique_children <- function(config_item_name) {
  data_syn_id <- config::get("data_syn_ids")[[config_item_name]]
  staging_syn_id <- config::get("staging_syn_ids")[[config_item_name]]

  children <- synGetChildren(data_syn_id)$asList() |>
    append(synGetChildren(staging_syn_id)$asList())

  # Get a list of files that appear in both Data and Staging so we can ignore
  # the Data version
  filenames <- sapply(children, "[[", "name")
  non_unique <- filenames[duplicated(filenames)] |> unique()

  final_list <- lapply(children, function(child) {
    if (child$name %in% non_unique) {
      child_info <- synGet(child$id, downloadFile = FALSE)
      if (child_info$parentId == data_syn_id) {
        message(str_glue("WARNING: `{child_info$name}` exists in both Data ",
                         "and Staging. Ignoring the Data version."))
        return(NULL)
      }
    }

    return(child)
  })

  final_list <- final_list[lengths(final_list) > 0]
  names(final_list) <- sapply(final_list, "[[", "name")

  return(final_list)
}


# Gets all harmonized metadata files from Synapse. Assumes that the latest
# version on Synapse is the desired file for each study.
#
# We have to check both "Staging/Metadata" and "Data/Metadata" for files, to
# account for cases where a file is in Staging and hasn't been curated and
# released yet.
get_all_metadata <- function() {
  meta_ids <- syn_get_unique_children("metadata")

  meta_list <- lapply(meta_ids, function(item) {
    syn_id <- paste0(item$id, ".", item$versionNumber)
    file_info <- synGet(syn_id, downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")

    data <- read.csv(file_info$path)

    return(list(
      study = unique(data$study),
      data = data,
      provenance = data.frame(
        id = file_info$id,
        versionNumber = file_info$versionNumber,
        name = file_info$name,
        study = unique(data$study)
      )
    ))
  })

  names(meta_list) <- sapply(meta_list, "[[", "study")
  return(meta_list)
}


# Download every count file for every study, read them in, and return them as a
# list. We have to check both "Staging/.../Raw Gene Counts" and
# "Data/.../Raw Gene Counts" for files, to account for cases where a file is in
# Staging and hasn't been curated and released yet.
#
# Arguments:
#   studies - a character vector of study names. Names must match the official
#     name of the study in the ADKP
#   meta_list - a list of metadata data frames, one per study
#   symbol_map - a data frame of Ensembl IDs and gene symbols for every gene in
#     the data set
#   count_type - one of "gene_counts", "gene_tpm", "transcript_counts", or
#     "transcript_tpm", corresponding to which count file to read in
#
# Returns:
#   a list of matrices, one per study
get_all_counts_files <- function(studies, meta_list, symbol_map,
                                 count_type = "gene_counts") {
  count_folders <- syn_get_unique_children("raw_gene_counts")
  names(count_folders) <- sapply(count_folders, "[[", "name")

  stopifnot(all(studies %in% names(count_folders)))

  counts_list <- lapply(studies, function(study_name) {
    print(study_name)
    folder <- count_folders[[study_name]]$id

    count_file_id <- synFindEntityId(str_glue("{study_name}.{count_type}.tsv"),
                                     parent = folder)

    if (is.null(count_file_id)) {
      warning(str_glue("Could not find a '{count_type}' file for {study_name}. ",
                       "This study will not be validated."))
      return(NULL)
    }

    counts_file <- synGet(count_file_id, downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

    counts <- read.delim(counts_file$path, header = TRUE, row.names = 1) |>
      # Get rid of transcript_id column
      select(-transcript_id.s.)

    # Replace specimenID column names with the unique_specimenID in the metadata
    id_map <- meta_list[[study_name]]
    rownames(id_map) <- id_map$R_safe_specimenID
    colnames(counts) <- id_map[colnames(counts), "unique_specimenID"]

    # Make sure the matrix has all genes and in the same order. Data sets that
    # were aligned to V1 of the genome will be missing human CLU, so this gene
    # is filled in with 0's.
    counts <- counts[symbol_map$ensembl_gene_id, ]
    rownames(counts) <- symbol_map$ensembl_gene_id
    counts[is.na(counts)] <- 0

    return(list(
      name = study_name,
      counts = counts,
      provenance = data.frame(
        id = counts_file$id,
        versionNumber = counts_file$versionNumber,
        name = counts_file$name,
        study = study_name
      )
    ))
  })

  # We can't add names() to the list because otherwise they get concatenated
  # to the column names when we cbind the list
  counts_list <- counts_list[lengths(counts_list) > 0]
  return(counts_list)
}


synapse_upload_de <- function(study_name, provenance_df, folder_syn_ids,
                              de_filenames, norm_filenames) {
  prov <- provenance_df[study_name, ]

  # There might be more than one DE or normalized file associated with the study
  for (filename in de_filenames) {
    syn_file <- File(filename, parent = folder_syn_ids$de_analysis)
    syn_file <- synStore(syn_file, forceVersion = FALSE,
                         set_annotations = FALSE,
                         used = prov$used, executed = prov$executed)
  }

  for (filename in norm_filenames) {
    syn_file <- File(filename, parent = folder_syn_ids$norm_counts) # TODO
    syn_file <- synStore(syn_file, forceVersion = FALSE,
                         set_annotations = FALSE,
                         used = prov$used, executed = prov$executed)
  }
}


# Check if a local file has an identical counterpart on Synapse by looking for
# the local file's Md5 sum
# Returns TRUE if a file with the same Md5 sum and filename exist on Synapse,
# FALSE otherwise.
syn_md5_check <- function(local_file) {
  md5 <- system(str_glue("md5sum {local_file}"), intern = TRUE) |>
    str_split("  ")
  md5 <- md5[[1]][1]

  syn_id <- synMd5Query(md5)

  return(length(syn_id) > 0 &&
           syn_id[[1]]$name == basename(local_file))
}


# Only uploads a file to Staging on Synapse if it doesn't exist in Data, or if
# it does exist but the version in Data doesn't have the same Md5 sum as the
# local file.
syn_safe_upload <- function(local_file, parent_id, used = NULL, executed = NULL) {
  # If a matching md5 sum is found on Synapse, return NULL
  if (syn_md5_check(local_file)) {
    return(NULL)
  }

  # Otherwise upload the file
  syn_file <- File(local_file, parent = parent_id)

  syn_file <- synStore(syn_file,
                       used = used,
                       executed = executed,
                       forceVersion = FALSE,
                       set_annotations = FALSE)

  print(str_glue("Stored {basename(local_file)}: ",
                 "{syn_file$id}.{syn_file$versionNumber}"))

  return(syn_file)
}


# Convert the list structure returned by synapserutils::walk() into a flat
# data frame for easier manipulation.
# `walk_list` should be the list returned by walk(...)$asList().
#
# The structure of each item in walk_list is as follows:
#   item[[1]] = list(1 = path of containing folder, 2 = synapse id of folder)
#   item[[2]] = list(empty) or list of sub-folder items (with the same fields as item[[1]])
#   item[[3]] = list of files, where each item is list(1 = file name, 2 = synapse id)
#
# Any sub-folders in item[[2]] also get their own entry in `walk_list`. We use
# this information to keep track of parent/child relationships.
syn_walk_to_df <- function(walk_list) {
  if (length(walk_list) == 0) {
    return(NULL)
  }

  walk_df <- lapply(walk_list, function(walk_item) {
    if (length(walk_item) < 3) {
      return(NULL)
    }

    sub_folders <- sapply(walk_item[[2]], "[[", 2)
    folder_df <- data.frame(name = basename(pluck(walk_item, 1, 1)),
                            id = pluck(walk_item, 1, 2),
                            path = pluck(walk_item, 1, 1),
                            children = paste(sub_folders, collapse = ", "),
                            type = "folder")
    file_df <- NULL

    if (length(walk_item[[3]]) != 0) {
      file_df <- lapply(walk_item[[3]], function(file_item) {
        data.frame(name = basename(file_item[[1]]),
                   id = file_item[[2]],
                   path = paste(folder_df$path, file_item[[1]], sep = "/"),
                   type = "file")
      }) |>
        list_rbind()

      folder_df$children <- paste(folder_df$children, file_df$id, collapse = ", ")
    }

    return(list_rbind(list(folder_df, file_df)))
  }) |>
    list_rbind()

  # Back-trace parent relationships
  walk_df$parent <- sapply(1:nrow(walk_df), function(N) {
    parent <- grepl(walk_df$id[N], walk_df$children)
    if (any(parent)) {
      walk_df$id[parent]
    } else {
      ""
    }
  })

  return(walk_df)
}
