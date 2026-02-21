get_all_metadata <- function(folder_synid) {
  meta_ids <- synGetChildren(folder_synid)$asList()

  meta_list <- lapply(meta_ids, function(item) {
    syn_id <- paste0(item$id, ".", item$versionNumber)
    file_info <- synGet(syn_id, downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")

    data <- read.csv(file_info$path)

    return(list(
      study = unique(data$study),
      data = data,
      provenance = data.frame(id = file_info$id,
                              versionNumber = file_info$versionNumber,
                              name = file_info$name)
    ))
  })

  names(meta_list) <- sapply(meta_list, "[[", "study")
  return(meta_list)
}


# Download every count file for every study, read them in, and return them as
# a list
#
# Arguments:
#   folder_synid - the Synapse ID of the raw counts folder that contains all
#     study count data
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
get_all_counts_files <- function(folder_synid, studies, meta_list, symbol_map,
                                 count_type = "gene_counts") {
  count_folders <- synGetChildren(folder_synid,
                                  includeTypes = list("folder"))$asList()
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

    return(list(name = study_name,
                provenance = paste0(counts_file$id, ".", counts_file$versionNumber),
                counts = counts))
  })

  # We can't add names() to the list because otherwise they get concatenated
  # to the column names when we cbind the list
  counts_list <- counts_list[lengths(counts_list) > 0]
  return(counts_list)
}
