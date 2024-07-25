# This script creates sample sheets from fastq files in the format needed for
# the Nextflow pipeline nf-core/rnaseq. This involves getting a list of fastq
# files from Synapse for each study and matching the files to the associated
# specimenID in the harmonized metadata. Each study formats their fastq
# filenames differently so there is a lot of study-specific handling of this.
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
# Note that although UCI_Trem2-R47H_NSS exists in the harmonized metadata file,
# there are some outstanding issues with the fastq filenames so they cannot be
# processed yet.
#
# The CSV file "Model_AD_SynID_list.csv" was created by hand and lists the
# Synapse IDs of all metadata files and the Synapse IDs of the folders containing
# fastq files. To add a new study to this process, a row needs to be added for
# that study in this CSV file, "01_Harmonize_Metadata.R" needs to be re-run to
# incorporate the metadata for that new study, and any special handling of fastq
# filenames needs to be added below.

library(synapser)
library(stringr)
library(dplyr)

synLogin()
tmp_dir <- file.path("data", "tmp")
samplesheet_dir <- file.path("data", "sample_sheets")

dir.create(tmp_dir, showWarnings = FALSE)
dir.create(samplesheet_dir, showWarnings = FALSE)

syn_ids <- read.csv(file.path("data", "Model_AD_SynID_list.csv"))

meta_file <- synGet("syn61850266", downloadLocation = tmp_dir)
metadata <- read.csv(file.path(tmp_dir, "Model_AD_merged_metadata.csv"))

# TEMPORARY: Skipping UCI_Trem2-R47H_NSS due to outstanding questions / issues
syn_ids <- subset(syn_ids, Study != "UCI_Trem2-R47H_NSS")


# Create one sample sheet per study --------------------------------------------

for (N in 1:nrow(syn_ids)) {
  row <- syn_ids[N,]
  print(row$Study)
  meta_filt <- subset(metadata, study_name == row$Study)

  fastq_folders <- str_split(row$Fastq_Folders, ";")[[1]]

  ## Get a list of all fastqs available ----------------------------------------

  all_fastqs <- lapply(fastq_folders, function(folder) {
    children <- synGetChildren(folder, includeTypes = list("file"))$asList()
    children <- do.call(rbind, children)
  })

  all_fastqs <- as.data.frame(do.call(rbind, all_fastqs))

  # Get which fastqs are read 1 and read 2
  # UCI_ABCA7 uses "1.fq" and "2.fq" instead of "R1.fastq" and "R2.fastq", and
  # the Jax studies use "R1_001.fastq.gz" and "R2_001.fastq.gz", so this regex
  # should capture all three cases.
  all_fastqs$read <- 1
  all_fastqs$read[grepl("2(_001)?\\.f[a|q]", all_fastqs$name)] <- 2


  ## Study-specific handling of filename matching ------------------------------

  if (row$Study %in% c("Jax.IU.Pitt_APOE4.Trem2.R47H", "Jax.IU.Pitt_5XFAD")) {
    # The first field in the filename is the individualID. "rh" needs to be
    # added at the end to get the specimenID.
    fastq_ids <- str_replace(all_fastqs$name, "_.*", "")
    all_fastqs$specimenID <- paste0(fastq_ids, "rh")

  } else if (row$Study %in% c("UCI_3xTg-AD", "UCI_Trem2_Cuprizone")) {
    # The individual ID is the field right before the sample number and read
    # number in the filename. For UCI_3xTg-AD, either "lh" or "rh" gets added
    # at the end to get the specimenID. For UCI_Trem2_Cuprizone, "b" is added
    # at the end.

    # Matches "_S#_R#.fastq.gz" where S# is something like S8, S12 and R# is R1 or R2
    fastq_ids <- str_replace(all_fastqs$name, "_S[0-9]+_R[1|2]\\.fastq\\.gz", "")
    fastq_ids <- str_replace(fastq_ids, ".*_", "")

    if (row$Study == "UCI_Trem2_Cuprizone") {
      all_fastqs$specimenID <- paste0(fastq_ids, "b")
    } else {
      # This is easier than figuring out whether to add "lh" or "rh"
      rownames(meta_filt) <- meta_filt$individualID
      all_fastqs$specimenID <- meta_filt[fastq_ids, "specimenID"]
    }

  } else if (row$Study == "UCI_5XFAD") {
    # There are two filename formats: For both, the individualID is in the first
    # field of the filename, but in one format it's followed by a C or H and in
    # the other it's preceded by "Hipp" or "Cortex". These values are used to
    # determine if "rc" or "rh" should be added to the individualID to get the
    # specimenID.
    fastq_ids <- str_split(all_fastqs$name, "_", simplify = TRUE)[, 1]
    cortex <- grepl("C", fastq_ids) # Catches both <id>C and Cortex<id>
    hipp <- grepl("H", fastq_ids) # Catches both <id>H and Hipp<id>

    # Change fastq name format to match specimen ID format
    fastq_ids[cortex] <- str_replace(fastq_ids[cortex], "Cortex", "")
    fastq_ids[cortex] <- str_replace(fastq_ids[cortex], "C", "")
    fastq_ids[cortex] <- paste0(fastq_ids[cortex], "rc")

    fastq_ids[hipp] <- str_replace(fastq_ids[hipp], "Hipp", "")
    fastq_ids[hipp] <- str_replace(fastq_ids[hipp], "H", "")
    fastq_ids[hipp] <- paste0(fastq_ids[hipp], "rh")

    all_fastqs$specimenID <- fastq_ids

  } else if (row$Study == "UCI_ABCA7") {
    # IndividualID is the field before the read number, which is either "1" or
    # "2" instead of "R1" or "R2". Files end in "fq.gz". "lh" is added to the
    # individualID to get the specimenID.
    fastq_ids <- str_replace(all_fastqs$name, "_[1|2]\\.fq\\.gz", "")
    fastq_ids <- str_replace(fastq_ids, ".*_", "")

    all_fastqs$specimenID <- paste0(fastq_ids, "lh")

  } else if (row$Study == "UCI_hAbeta_KI") {
    # IndividualID is the first field in the filename, including "-" characters.
    # Some IDs need to be altered to match what's in the metadata. "rh" is added
    # to the end of the individualID to get specimenID.
    fastq_ids <- str_replace(all_fastqs$name, "_.*", "")
    fastq_ids <- str_replace(fastq_ids, "C57-", "C57")
    fastq_ids <- str_replace(fastq_ids, "41-2-4", "41-24")

    all_fastqs$specimenID <- paste0(fastq_ids, "rh")

  } else if (row$Study == "UCI_PrimaryScreen") {
    # There are 3 different naming formats for these files:
    #   1. The individualID is first in the string
    #   2. The individualID is just before "_S#_R1/2.fastq.gz"
    #   3. The individualID is just before "_R1/2.fastq.gz"
    # After extracting the individualID from each filename, "rh" is added to
    # the end to get specimenIDs.

    # Case 1 handling -- files that start with a 4-digit number
    case1 <- grepl("^[0-9]{4}", all_fastqs$name)
    case1_ids <- str_replace(all_fastqs$name[case1], "_.*", "")

    # Case 2 handling and case 3 handling -- all non-case1 files.
    # Removing the read number and then the sample number makes case 2 and 3
    # handling identical.
    case2_3 <- !grepl("^[0-9]{4}", all_fastqs$name)
    case2_3_ids <- str_replace(all_fastqs$name[case2_3], "_R[1|2].*", "")
    case2_3_ids <- str_replace(case2_3_ids, "_S[0-9]+$", "")
    case2_3_ids <- str_replace(case2_3_ids, ".*_", "")

    all_fastqs$specimenID <- ""
    all_fastqs$specimenID[case1] <- case1_ids
    all_fastqs$specimenID[case2_3] <- case2_3_ids
    all_fastqs$specimenID <- paste0(all_fastqs$specimenID, "rh")

  } else if (row$Study == "UCI_Trem2-R47H_NSS") {
    # Matches "_S#_R#.fastq.gz" where S# is something like S8, S12 and R# is R1 or R2
    fastq_ids <- str_replace(all_fastqs$name, "_S[0-9]+_R[1|2]\\.fastq\\.gz", "")
    fastq_ids <- str_replace(fastq_ids, ".*_", "")
    # TODO waiting on answers about duplicate/mis-named fastq files
  }

  ## Format sample sheet for NextFlow ------------------------------------------

  lines <- all_fastqs %>%
            group_by(specimenID) %>%
            summarize(fastq_1 = paste0("syn://", id[read == 1]),
                      fastq_2 = paste0("syn://", id[read == 2]),
                      strandedness = "auto",
                      .groups = "drop") %>%
            rename(sample = specimenID)

  write.csv(lines,
            file.path(samplesheet_dir, paste0(row$Study, "_samplesheet_synapse.csv")),
            row.names = FALSE, quote = FALSE)


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
