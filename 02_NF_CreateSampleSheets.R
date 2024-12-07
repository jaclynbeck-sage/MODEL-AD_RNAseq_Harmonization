# This script creates sample sheets from fastq files in the format needed for
# the Nextflow pipeline nf-core/rnaseq. This involves getting a list of fastq
# files from Synapse for each study and matching the files to the associated
# specimenID in the harmonized metadata. Most fastq files are annotated with
# the specimenID but some studies need extra handling for mismatches.
#
# The following studies are currently accounted for in this script:
#   Jax.IU.Pitt_5XFAD
#   Jax.IU.Pitt_APOE4.Trem2.R47H
#   UCI_3xTg-AD
#   UCI_5XFAD
#   UCI_ABCA7
#   UCI_hAbeta_KI
#   UCI_PrimaryScreen
#   UCI_Trem2_Cuprizone
#
# The CSV file "Model_AD_SynID_list.csv" was created by hand and lists the
# Synapse IDs of all metadata files and the Synapse IDs of the folders containing
# fastq files. To add a new study to this process, a row needs to be added for
# that study in this CSV file, "01_Harmonize_Metadata.R" needs to be re-run to
# incorporate the metadata for that new study, and any special handling of fastq
# filenames needs to be added below.
#
# After running this script, sample sheets were uploaded to Sage's Nextflow
# Tower environment, and the nf-synapse pipeline was run to download all of
# the fastq files and convert the sample sheets to point to their locations on
# the S3 bucket of the environment. Then, the nf-core/rnaseq pipeline was run
# with the updated sample sheets to align the data.

library(synapser)
library(stringr)
library(dplyr)

# Synapse IDs used in this script that are not in Model_AD_SynID_list.csv
syn_metadata_file_id <- "syn61850266"
syn_ref_fasta_id <- "syn62035247"
syn_ref_gtf_id <- "syn62035250"
syn_jax_5xfad_fastq_key_id <- "syn34733116"
syn_samplesheet_folder_id <- "syn62147112"
syn_portal_query_id <- "syn11346063"

synLogin(silent = TRUE)
tmp_dir <- file.path("data", "tmp")
samplesheet_dir <- file.path("data", "sample_sheets")
provenance_dir <- file.path("data", "provenance_manifests")

dir.create(tmp_dir, showWarnings = FALSE)
dir.create(samplesheet_dir, showWarnings = FALSE)
dir.create(provenance_dir, showWarnings = FALSE)

syn_ids <- read.csv(file.path("data", "Model_AD_SynID_list.csv"))

meta_file <- synGet(syn_metadata_file_id, downloadLocation = tmp_dir,
                    ifcollision = "overwrite.local")
metadata <- read.csv(meta_file$path)

ref_fasta <- synGet(syn_ref_fasta_id, downloadFile = FALSE)
ref_gtf <- synGet(syn_ref_gtf_id, downloadFile = FALSE)

meta_provenance <- rbind(
  c(meta_file$id, meta_file$versionNumber, meta_file$name),
  c(ref_fasta$id, ref_fasta$versionNumber, ref_fasta$name),
  c(ref_gtf$id, ref_gtf$versionNumber, ref_gtf$name)
)
colnames(meta_provenance) <- c("id", "currentVersion", "name")


# Create one sample sheet per study --------------------------------------------

for (N in 1:nrow(syn_ids)) {
  row <- syn_ids[N,]
  print(row$Study)
  meta_filt <- subset(metadata, study_name == row$Study)

  ## Get a list of all fastqs available ----------------------------------------

  # Query the portal -- this query will work correctly for all studies except
  # UCI_ABCA7
  query <- paste0("SELECT * FROM ", syn_portal_query_id, " WHERE ",
                  "( ( \"assay\" HAS ( 'long-read rnaSeq', 'rnaSeq' ) ) AND ",
                  "( \"fileFormat\" = 'fastq' ) AND ",
                  "( \"isMultiSpecimen\" = 'false' OR \"isMultiSpecimen\" IS NULL ) AND ",
                  "( \"study\" HAS ( '", row$Study, "' ) ) )")

  # UCI_ABCA7 fastq files are not all properly annotated
  if (row$Study == "UCI_ABCA7") {
    query <- paste0("SELECT * FROM ", syn_portal_query_id, " WHERE ",
                    "( ( \"name\" LIKE '%fq.gz%' ) ) AND ",
                    "( ( ( \"study\" HAS ( '", row$Study, "' ) ) ) )")
  }

  result <- synTableQuery(query, includeRowIdAndRowVersion = FALSE)

  all_fastqs <- read.csv(result$filepath) %>%
    select(id, name, specimenID, currentVersion)

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

  # The annotated specimenIDs don't match the updated metadata files for the
  # Jax studies so we map between the two using the fastq files key provided
  # with the updated metadata
  if (row$Study %in% c("Jax.IU.Pitt_5XFAD", "Jax.IU.Pitt_APOE4.Trem2.R47H")) {
    fastq_key_id <- ifelse(row$Study == "Jax.IU.Pitt_5XFAD",
                           syn_jax_5xfad_fastq_key_id,
                           row$Metadata_Assay)

    fastq_key_file <- synGet(fastq_key_id, downloadLocation = tmp_dir,
                             ifcollision = "overwrite.local")
    fastq_key_df <- read.csv(fastq_key_file$path) %>%
      select(sampleName, fastq_1, fastq_2) %>%
      tidyr::pivot_longer(cols = c(fastq_1, fastq_2),
                          names_to = "fastq_type",
                          values_to = "name") %>%
      mutate(name = str_replace(name, "_S.*_L00.", "")) %>%
      distinct() %>%
      dplyr::rename(specimenID = sampleName)

    all_fastqs <- all_fastqs %>%
      select(-specimenID) %>%
      merge(fastq_key_df)

  } else if (row$Study == "UCI_5XFAD") {
    # Specimen IDs in the annotation are formatted with a numerical ID followed
    # by "C_RNAseq" or "H_RNAseq", e.g. "305C_RNAseq". However in the metadata
    # files the IDs are formatted as the numerical ID followed by "rc" or "rh",
    # so this converts to the right format.
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "C_.*", "rc")
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "H_.*", "rh")

  } else if (row$Study == "UCI_ABCA7") {
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

  } else if (row$Study == "UCI_hAbeta_KI") {
    # Remove some special characters (commas and parentheses) to match what's in
    # the metadata
    all_fastqs$specimenID <- str_replace_all(all_fastqs$specimenID,
                                             ",|\\(|\\)", "")

  } else if (row$Study == "UCI_Trem2_Cuprizone") {
    # Some specimenIDs have an extra "w" in them
    all_fastqs$specimenID <- str_replace(all_fastqs$specimenID, "w", "")

  } else if (row$Study == "UCI_Trem2-R47H_NSS") {
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


  ## Format sample sheet for NextFlow ------------------------------------------

  lines <- all_fastqs %>%
            group_by(specimenID) %>%
            summarize(fastq_1 = paste0("syn://", id[read == 1]),
                      fastq_2 = paste0("syn://", id[read == 2]),
                      strandedness = "auto",
                      .groups = "drop") %>%
            dplyr::rename(sample = specimenID)

  samplesheet_filename <- file.path(samplesheet_dir,
                                    paste0(row$Study, "_samplesheet_rnaseq.csv"))
  write.csv(lines, samplesheet_filename, row.names = FALSE, quote = FALSE)


  ## Create provenance manifest for Step 03 ------------------------------------

  provenance <- select(all_fastqs, id, currentVersion, name)
  provenance <- rbind(meta_provenance, provenance)

  write.csv(provenance,
            file.path(provenance_dir, paste0(row$Study, "_provenance.csv")),
            row.names = FALSE, quote = FALSE)


  ## Upload sample sheet to Synapse --------------------------------------------

  # Remove universal reference genome IDs as they are not used for sample sheets,
  # only in step 03.
  samplesheet_provenance <- subset(provenance,
                                   !grepl("(fa|gtf)\\.gz", provenance$name))
  used_ids <- paste(samplesheet_provenance$id,
                    samplesheet_provenance$currentVersion,
                    sep = ".")

  github_link <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/02_NF_CreateSampleSheets.R"
  syn_file <- File(samplesheet_filename, parent = syn_samplesheet_folder_id)
  synStore(syn_file,
           used = used_ids,
           executed = github_link,
           forceVersion = FALSE)


  ## Print warnings for missing fastqs or extra fastqs -------------------------

  if (!all(all_fastqs$specimenID %in% meta_filt$specimenID)) {
    extra <- setdiff(all_fastqs$specimenID, meta_filt$specimenID)
    msg <- paste("WARNING:", row$Study, "has", length(extra), "fastq files for",
                 "specimens that do not exist in the metadata:",
                 paste(sort(extra), collapse = ", "))
    message(msg)
  }

  if (!all(meta_filt$specimenID %in% all_fastqs$specimenID)) {
    missing <- setdiff(meta_filt$specimenID, all_fastqs$specimenID)
    msg <- paste("WARNING:", row$Study, "is missing fastq files for",
                 length(missing), "specimens:",
                 paste(sort(missing), collapse = ", "))
    message(msg)
  }
}
