# This script updates the RSEM output file names on Synapse and adds provenance.
# The counts files for every study are by default named
# "rsem.merged.<count_type>.tsv", which could cause confusion or accidental
# over-writes when downloading the data. Here we replace "rsem.merged" in each
# file name with the name of the study, and add provenance pointing to the
# metadata, fastq files, genome reference, and NextFlow pipeline input files.
#
# This script is a little fragile because of string matching and assumes:
#   1. Step 03 has been run in the same directory so that provenance files exist
#      in the "output/provenance_manifests" directory.
#   2. The NextFlow config and samplesheet files contain the study name in the
#      filename exactly as specified in the Model_AD_SynID_list.csv file
#   3. The RSEM counts files all exist in a folder under syn51132850, and that
#      folder is named with the study exactly as specified in the syn ID list
#      file.

library(synapser)
library(synapserutils)
library(stringr)

# TODO eventually maybe most of this string matching can be done with annotations

folder_syn_ids <- config::get("folder_syn_ids", config = "default")

synLogin(silent = TRUE)

provenance_dir <- file.path("output", "provenance_manifests")
syn_id_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"),
                        comment.char = "#")

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

# Get NextFlow configuration files and samplesheet files
nf_config_files <- synGetChildren(folder_syn_ids$nf_configs)$asList()
nf_samplesheet_files <- synGetChildren(folder_syn_ids$nf_samplesheets)$asList()

nf_config_names <- sapply(nf_config_files, "[[", "name")
nf_samplesheet_names <- sapply(nf_samplesheet_files, "[[", "name")


# Loop through each study folder -----------------------------------------------

for (folder in study_folders) {
  study <- folder$name
  print(study)

  # Get all files in the study folder. The only files in there should be RSEM
  # output files. This will skip folders (like the quality control folder).
  syn_files <- synGetChildren(folder$id,
                              includeTypes = list("file"))
  syn_files <- syn_files$asList()


  ## Rename each file and add provenance ---------------------------------------

  for (sf in syn_files) {
    # If we haven't already renamed this file, rename it. The files are
    # originally called "rsem.merged.gene_counts.tsv", "rsem.merged.gene_tpm.tsv",
    # "rsem.merged.transcript_counts.tsv", or "rsem.merged.transcript_tpm.tsv".
    # Here we replace "rsem.merged" with the study name.
    if (!grepl(study, sf$name)) {
      new_name <- str_replace(sf$name, "rsem\\.merged", study)
      print(str_glue("Renaming {sf$name} to {new_name}"))

      changeFileMetaData(sf$id,
                         downloadAs = new_name,
                         name = new_name,
                         forceVersion = FALSE)
    }


    ### Check that necessary files exist ---------------------------------------

    prov_file <- file.path(provenance_dir,
                            paste0(study, "_provenance.csv"))

    # Check for the existance of the provinance manifest, NextFlow config,
    # and samplesheet that match this study. The provenance manifest was
    # created in step 02.
    if (!file.exists(prov_file)) {
      msg <- paste0("WARNING: cannot find provenance manifest for ",
                    study, ". Provenance will not be set.")
      message(msg)
      next
    }

    config_match <- which(grepl(paste0("rnaseq_", study),
                                nf_config_names))
    samplesheet_match <- which(grepl(paste0(study, ".*_rnaseq"),
                                     nf_samplesheet_names))

    if (length(config_match) == 0 || length(samplesheet_match) == 0) {
      msg <- paste0("WARNING: either the NextFlow configuration file or the ",
                    "sample sheet file are missing from Synapse for ", study,
                    ". Provenance will not be set.")
      message(msg)
      next
    }


    ### Set provenance on the file ---------------------------------------------

    provenance <- read.csv(prov_file)
    nf_config <- nf_config_files[[config_match]]
    nf_samplesheet <- nf_samplesheet_files[[samplesheet_match]]

    # the "versionNumber" column is called "currentVersion" in the CSV file due
    # to it being named that way in the synTableQuery in Step 02. The column
    # will stay named "currentVersion" despite using $versionNumber in the two
    # new rows.
    provenance <- rbind(
      c(nf_config$id, nf_config$versionNumber, nf_config$name),
      c(nf_samplesheet$id, nf_samplesheet$versionNumber, nf_samplesheet$name),
      provenance
    )

    act <- Activity(used = paste(provenance$id,
                                 provenance$currentVersion,
                                 sep = "."))
    synSetProvenance(sf$id, act)
  }
}
