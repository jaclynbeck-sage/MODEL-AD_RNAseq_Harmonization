# This script creates sample sheets from fastq files in the format needed for
# the Nextflow pipeline nf-core/rnaseq. This involves getting a list of fastq
# files from Synapse for each study and matching the files to the associated
# specimenID in the harmonized metadata. Most fastq files are annotated with the
# specimenID but some studies need extra handling for mismatches.
#
# The following studies are currently accounted for in this script:
#   Jax.IU.Pitt_5XFAD
#   Jax.IU.Pitt_APOE4.Trem2.R47H
#   Jax.IU.Pitt_LOAD2.PrimaryScreen
#   UCI_3xTg-AD
#   UCI_5XFAD
#   UCI_ABCA7
#   UCI_Bin1K358R
#   UCI_Clu-h2kbKI
#   UCI_hAbeta_KI
#   UCI_PrimaryScreen
#   UCI_Trem2_Cuprizone
#
# To add a new study, make sure steps 01 and 02 have been run for the new study,
# and ensure that the study name is in the `studies` field of the `config.yml`
# file. Add any special handling of fastq filenames for the new study below.
#
# After running this script, sample sheets were uploaded to Sage's Nextflow
# Tower environment, and the nf-synapse pipeline was run to download all of the
# fastq files and convert the sample sheets to point to their locations on the
# S3 bucket of the environment. Then, the nf-core/rnaseq pipeline was run with
# the updated sample sheets to align the data.

library(synapser)
library(stringr)
library(dplyr)

source("util_functions.R")

file_syn_ids <- config::get("file_syn_ids")
staging_syn_ids <- config::get("staging_syn_ids")
syn_portal_query_id <- config::get("adkp_query_id")
studies <- config::get("studies")

synLogin(silent = TRUE)
tmp_dir <- file.path("output", "tmp")
samplesheet_dir <- file.path("output", "sample_sheets")
provenance_dir <- file.path("output", "provenance_manifests")

dir.create(tmp_dir, showWarnings = FALSE)
dir.create(samplesheet_dir, showWarnings = FALSE)
dir.create(provenance_dir, showWarnings = FALSE)

meta_list <- get_all_metadata()

ref_fasta <- synGet(file_syn_ids$ref_fasta, downloadFile = FALSE)
ref_gtf <- synGet(file_syn_ids$ref_gtf, downloadFile = FALSE)

meta_provenance <- rbind(
  c(ref_fasta$id, ref_fasta$versionNumber, ref_fasta$name),
  c(ref_gtf$id, ref_gtf$versionNumber, ref_gtf$name)
)
colnames(meta_provenance) <- c("id", "versionNumber", "name")


# Create one sample sheet per study --------------------------------------------

for (study in studies) {
  meta <- meta_list[[study]]$data

  ## Get a list of all fastqs available ----------------------------------------

  # Query the portal -- this query will work correctly for all studies
  query <- str_glue(
    "SELECT * FROM {syn_portal_query_id} WHERE ",
    "( ( \"assay\" HAS ( 'long-read rnaSeq', 'rnaSeq' ) ) AND ",
    "( \"fileFormat\" = 'fastq' ) AND ",
    "( \"isMultiSpecimen\" = 'false' OR \"isMultiSpecimen\" IS NULL ) AND ",
    "( \"study\" HAS ( '{study}' ) ) )"
  )

  result <- synTableQuery(query, includeRowIdAndRowVersion = FALSE)

  all_fastqs <- read.csv(result$filepath) |>
    select(id, name, specimenID, currentVersion) |>
    dplyr::rename(versionNumber = currentVersion)

  # Get which fastqs are read 1 and read 2
  # UCI_ABCA7 uses "1.fq" and "2.fq" instead of "R1.fastq" and "R2.fastq", and
  # the Jax studies use "R1_001.fastq.gz" and "R2_001.fastq.gz", so this regex
  # should capture all three cases. It recognizes "2", which may or may not be
  # followed by "_001", followed by ".fa" or ".fq".
  all_fastqs$read <- 1
  all_fastqs$read[grepl("2(_001)?\\.f[a|q]", all_fastqs$name)] <- 2

  # Check that we have an equal number of read 1 and read 2
  n_reads <- table(all_fastqs$read)
  stopifnot(n_reads[1] == n_reads[2])


  ## Study-specific handling of ID formatting issues ---------------------------

  # Temporary: The annotated specimenIDs don't match the updated metadata files
  # for the Jax studies, but the specimenID is in the filename
  if (study == "Jax.IU.Pitt_5XFAD") {
    tmp_ids <- str_split(all_fastqs$name, "_", simplify = TRUE)[, c(2, 3)]
    all_fastqs$specimenID <- paste0(tmp_ids[, 1], "_", tmp_ids[, 2])

  } else if (study == "UCI_5XFAD") {
    # Specimen IDs in the annotation are formatted with a numerical ID followed
    # by "C_RNAseq" or "H_RNAseq", e.g. "305C_RNAseq". However in the metadata
    # files the IDs are formatted as the numerical ID followed by "rc" or "rh",
    # so this converts to the right format.
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "C_.*", "rc")
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "H_.*", "rh")

  } else if (study == "UCI_ABCA7") {
    # There are 4 fastq files that don't have an annotated specimen ID, but
    # the filename contains an ID that exists in the assay metadata. For these
    # 4 files, we extract the ID from the name and add "lh" to it to get the
    # specimen ID.
    missing_inds <- which(is.na(all_fastqs$specimenID) |
                            nchar(all_fastqs$specimenID) == 0)
    fastq_names <- all_fastqs$name[missing_inds]

    # IndividualID is the field before the read number, which is either "1" or
    # "2" instead of "R1" or "R2". Files end in "fq.gz". "lh" is added to the
    # individualID to get the specimenID.
    fastq_ids <- str_replace(fastq_names, "_[1|2]\\.fq\\.gz", "")
    fastq_ids <- str_replace(fastq_ids, ".*_", "")

    all_fastqs$specimenID[missing_inds] <- paste0(fastq_ids, "lh")

  } else if (study == "UCI_Clu-h2kbKI") {
    # Fastqs for specimenID "12680lc" are in a folder marked "Deprecated" and
    # should be excluded. There are 8 more fastq files in that folder for 4
    # specimenIDs that don't exist in the metadata, so these are excluded as
    # well.
    all_fastqs <- subset(all_fastqs, specimenID != "12680lc" &
                           specimenID %in% meta$specimenID)

  } else if (study == "UCI_Trem2_Cuprizone") {
    # Some specimenIDs have an extra "w" in them
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "w", "")

  } else if (study == "UCI_Trem2-R47H_NSS") {
    # Some samples were re-run so there are 4 associated fastq files instead of
    # 2. In this case, we want to use only the re-sequenced files, which will
    # have "RB1" or "RB2" in the name.
    id_tbl <- table(all_fastqs$specimenID)
    duplicate_ids <- names(id_tbl)[id_tbl > 2]

    all_fastqs$is_duplicate <- all_fastqs$specimenID %in% duplicate_ids
    all_fastqs$is_reseq <- grepl("_RB[1|2]_", all_fastqs$name)

    keep <- (!all_fastqs$is_duplicate) | all_fastqs$is_reseq
    all_fastqs <- all_fastqs[keep, ]
  }

  # Check for fastq files for specimens not in the metadata, and check for
  # duplicate fastq files. Each specimenID should only have 2 files.
  stopifnot(all_fastqs$specimenID %in% meta$specimenID)
  stopifnot(all(table(all_fastqs$specimenID) == 2))


  ## Format sample sheet for NextFlow ------------------------------------------

  lines <- all_fastqs |>
            group_by(specimenID) |>
            summarize(fastq_1 = paste0("syn://", id[read == 1]),
                      fastq_2 = paste0("syn://", id[read == 2]),
                      strandedness = "auto",
                      .groups = "drop") |>
            dplyr::rename(sample = specimenID)

  samplesheet_filename <- file.path(samplesheet_dir,
                                    paste0(study, "_samplesheet_rnaseq.csv"))
  write.csv(lines, samplesheet_filename, row.names = FALSE, quote = FALSE)


  ## Create provenance manifest for Step 03 ------------------------------------

  provenance <- select(all_fastqs, id, versionNumber, name)
  provenance <- do.call(rbind,
                        list(select(meta_list[[study]]$provenance, -study),
                             meta_provenance,
                             provenance))

  write.csv(provenance,
            file.path(provenance_dir, paste0(study, "_provenance.csv")),
            row.names = FALSE, quote = FALSE)


  ## Upload sample sheet to Synapse --------------------------------------------

  # Remove universal reference genome IDs as they are not used for sample sheets,
  # only in step 03.
  samplesheet_provenance <- subset(provenance,
                                   !grepl("(fa|gtf)\\.gz", provenance$name))
  used_ids <- paste(samplesheet_provenance$id,
                    samplesheet_provenance$versionNumber,
                    sep = ".")

  github_link <- paste0(config::get("github_repo_url"),
                        "/blob/main/03_NF_Create_FASTQ_Sample_Sheets.R")

  syn_safe_upload(samplesheet_filename,
                  parent_id = staging_syn_ids$nf_samplesheets,
                  used = used_ids,
                  executed = github_link)


  ## Print warnings for missing fastqs or extra fastqs -------------------------

  if (!all(all_fastqs$specimenID %in% meta$specimenID)) {
    extra <- setdiff(all_fastqs$specimenID, meta$specimenID)
    msg <- str_glue(
      "WARNING: {study} has {length(extra)} fastq files for specimen(s) ",
      "that do not exist in the metadata: {paste(sort(extra), collapse = ", ")}"
    )
    message(msg)
  }

  if (!all(meta$specimenID %in% all_fastqs$specimenID)) {
    missing <- setdiff(meta$specimenID, all_fastqs$specimenID)
    msg <- str_glue(
      "WARNING: {study} is missing fastq files for {length(missing)} ",
      "specimen(s): {paste(sort(missing), collapse = ", ")}"
    )
    message(msg)
  }
}
