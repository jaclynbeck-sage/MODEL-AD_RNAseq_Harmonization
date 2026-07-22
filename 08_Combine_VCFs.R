# This script loads the VCF files (output by nf-core/sarek) for each study and
# concatenates the data into a single CSV file per study. The CSVs are uploaded
# to the Genotype Validation folder in Synapse.
#
# Locating, downloading, and reading all of the VCF files can take several
# minutes per study, which is why this is a separate step from actual genotype
# validation.

library(synapser)
library(synapserutils)
library(stringr)
library(dplyr)

source("util_functions.R")

# Set up -----------------------------------------------------------------------

# Set to TRUE if only VCFs in Staging/ should be combined, i.e. all released
# data in Data/ is ignored. Set to FALSE to combine VCFs from all data in both
# Staging/ and Data/.
staging_only <- TRUE

file_syn_ids <- config::get("file_syn_ids")
staging_syn_ids <- config::get("staging_syn_ids")
studies <- config::get("studies")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url"),
                 "/blob/main/08_Combine_VCFs.R")
tmp_dir <- file.path("output", "tmp")

vcf_dir <- file.path("output", "combined_vcfs")
dir.create(vcf_dir, showWarnings = FALSE)


# Find vcf files on Synapse ----------------------------------------------------

# This data has a number of sub-folders and sub-sub-folders, so we use the walk
# function from synapserutils to get all the files. We use syn_get_unique_children
# to get all the vcf folders in both Data/ and Staging/ if `staging_only` is
# FALSE.

if (staging_only) {
  vcf_folders <- synGetChildren(staging_syn_ids$genotype_validation)$asList()
} else {
  vcf_folders <- syn_get_unique_children("genotype_validation")
}

names(vcf_folders) <- sapply(vcf_folders, "[[", "name")

vcf_folders <- vcf_folders[names(vcf_folders) %in% studies]

cat("Finding VCF files on Synapse...\n")
study_vcfs <- lapply(vcf_folders, function(folder) {
  print(folder$name)
  folder_structure <- synapserutils::walk(folder$id,
                                          includeTypes = list("file"))$asList()

  samples_df <- syn_walk_to_df(folder_structure) |>
    # Only keep files that end in .gz, not .gz.tbi
    subset(type == "file" & grepl("gz$", name)) |>
    mutate(study = folder$name,
           # Specimen ID is the name of the containing folder
           specimenID = basename(dirname(path))) |>
    select(name, id, specimenID, study) |>
    dplyr::rename(filename = name, syn_id = id)

  return(samples_df)
})

study_vcfs <- study_vcfs[lengths(study_vcfs) > 0]


# Download and parse vcf files -------------------------------------------------

cat("Reading VCF files...\n")
for (study_df in study_vcfs) {
  print(unique(study_df$study))

  mutation_df <- apply(study_df, 1, function(row) {
    vcf_file <- synGet(
      row[["syn_id"]],
      downloadLocation = file.path(tmp_dir, unique(study_df$study)),
      ifcollision = "overwrite.local"
    )

    vcf <- NULL

    # If the sample didn't have any detected mutations, the vcf file will
    # contain only comments and no data, which will cause the read.delim
    # function to throw an error.
    try({
      vcf <- read.delim(gzfile(vcf_file$path), header = FALSE, comment.char = "#") |>
        dplyr::rename(chromosome = V1, position = V2, ref = V4, alt = V5,
                      quality = V6) |>
        select(chromosome, position, ref, alt, quality) |>
        mutate(study = row[["study"]],
               specimenID = row[["specimenID"]])
    }, silent = TRUE)

    # If no mutations found, this will be NULL
    return(vcf)
  })

  mutation_df <- do.call(rbind, mutation_df)

  output_file <- file.path(vcf_dir,
                           str_glue("{unique(study_df$study)}_combined_vcf.csv"))
  write.csv(mutation_df, output_file, quote = FALSE, row.names = FALSE)

  # Upload to Synapse Staging
  syn_file <- syn_safe_upload(output_file, staging_syn_ids$genotype_validation,
                              used = study_df$syn_id, executed = github)
}
