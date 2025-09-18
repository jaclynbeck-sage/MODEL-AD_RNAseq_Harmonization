# This script creates sample sheets from BAM files in the format needed for
# the Nextflow pipeline nf-core/sarek. This involves getting a list of BAM
# files from Synapse for each study and matching the files to the associated
# specimenID in the harmonized metadata.
#
# This script assumes that there is a folder under syn63856101 whose name
# exactly matches each study name, and that matching folder has BAM files
# inside that are named in the format <specimenID>.merged.markdup.sorted.bam[.bai]
#
# After running this script, sample sheets were uploaded to Sage's Nextflow
# Tower environment, and the nf-synapse pipeline was run to download all of
# the BAM files and convert the sample sheets to point to their locations on
# the S3 bucket of the environment. Then, the nf-core/sarek pipeline was run
# with the updated sample sheets to call genotypes for each sample.

library(synapser)
library(stringr)
library(dplyr)

file_syn_ids <- config::get("file_syn_ids", config = "default")
folder_syn_ids <- config::get("folder_syn_ids", config = "default")

synLogin(silent = TRUE)
tmp_dir <- file.path("output", "tmp")
samplesheet_dir <- file.path("output", "sample_sheets")
provenance_dir <- file.path("output", "provenance_manifests")

dir.create(tmp_dir, showWarnings = FALSE)
dir.create(samplesheet_dir, showWarnings = FALSE)
dir.create(provenance_dir, showWarnings = FALSE)

meta_file <- synGet(file_syn_ids$merged_metadata, downloadLocation = tmp_dir,
                    ifcollision = "overwrite.local")
metadata <- read.csv(meta_file$path)

meta_provenance <- c(meta_file$id, meta_file$versionNumber, meta_file$name)
names(meta_provenance) <- c("id", "versionNumber", "name")

bam_folders <- synGetChildren(folder_syn_ids$bams)$asList()
names(bam_folders) <- sapply(bam_folders, "[[", "name")


# Create one sample sheet per study --------------------------------------------

for (study in unique(metadata$study_name)) {
  print(study)
  meta_filt <- subset(metadata, study_name == study) %>%
    select(individualID, specimenID)

  ## Get a list of all bam files -----------------------------------------------
  bam_files <- synGetChildren(bam_folders[[study]])$asList()

  if (length(bam_files) == 0) {
    message(str_glue("WARNING: {study} BAM folder is empty. Skipping..."))
    next
  }

  bam_files <- as.data.frame(do.call(rbind, bam_files)) %>%
    select(name, id, versionNumber)

  # specimenID is everything before ".markdup"
  bam_files$specimenID <- str_replace(bam_files$name, "\\.markdup.*", "")

  # file_type is the last 3 characters
  bam_files$file_type <- str_sub(bam_files$name, start = -3)

  # This will include BAM files with specimenIDs that don't exist in the
  # metadata, but the individualID will be NA
  bam_files <- merge(bam_files, meta_filt, by = "specimenID", all.x = TRUE) %>%
    subset(file_type %in% c("bam", "bai"))

  # Check that we have an equal number of bams and bais
  n_bams <- table(bam_files$file_type)
  stopifnot(n_bams[1] == n_bams[2])


  ## Format sample sheet for NextFlow ------------------------------------------

  lines <- bam_files %>%
    group_by(individualID, specimenID) %>%
    summarize(bam = paste0("syn://", id[file_type == "bam"]),
              bai = paste0("syn://", id[file_type == "bai"]),
              .groups = "drop") %>%
    dplyr::rename(patient = individualID, sample = specimenID)

  samplesheet_filename <- file.path(samplesheet_dir,
                                    paste0(study, "_samplesheet_variant_calling.csv"))
  write.csv(lines, samplesheet_filename, row.names = FALSE, quote = FALSE)


  ## Create provenance manifest for Step 03 ------------------------------------

  provenance <- select(bam_files, id, versionNumber, name) %>%
    mutate(across(id:name, unlist))

  provenance <- rbind(meta_provenance, provenance)

  write.csv(provenance,
            file.path(provenance_dir,
                      paste0(study, "_provenance_variant_calling.csv")),
            row.names = FALSE, quote = FALSE)


  ## Upload sample sheet to Synapse --------------------------------------------

  used_ids <- paste(provenance$id,
                    provenance$versionNumber,
                    sep = ".")

  github_link <- paste0(config::get("github_repo_url", config = "default"),
                        "/blob/main/05_NF_Create_BAM_Sample_Sheets.R")
  syn_file <- File(samplesheet_filename, parent = folder_syn_ids$nf_samplesheets)
  synStore(syn_file,
           used = used_ids,
           executed = github_link,
           forceVersion = FALSE)


  ## Print warnings for missing fastqs or extra fastqs -------------------------

  if (!all(bam_files$specimenID %in% meta_filt$specimenID)) {
    extra <- setdiff(bam_files$specimenID, meta_filt$specimenID)
    msg <- str_glue(
      "WARNING: {study} has BAM files for {length(extra)} specimen(s) that do ",
      "not exist in the metadata: {paste(sort(extra), collapse = ", ")}"
    )
    message(msg)
  }

  if (!all(meta_filt$specimenID %in% bam_files$specimenID)) {
    missing <- setdiff(meta_filt$specimenID, bam_files$specimenID)
    msg <- str_glue(
      "WARNING: {study} is missing BAM files for {length(missing)} ",
      "specimen(s): {paste(sort(missing), collapse = ", ")}"
    )
    message(msg)
  }
}
