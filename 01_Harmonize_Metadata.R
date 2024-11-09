# This script pulls in all metadata files necessary for RNA seq analysis from
# all studies included in the analysis, cleans up a few fields, and combines the
# metadata into one data frame containing all study information. Generally 3
# metadata files are needed: an individual metadata file, a biospecimen metadata
# file, and an RNA seq assay metadata file, which are merged together so that
# the final data frame only contains rows for specimens that have RNA seq data.
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
#   UCI_Trem2-R47H_NSS
#
# The CSV file "Model_AD_SynID_list.csv" was created by hand and lists the
# Synapse IDs of all metadata files and the Synapse IDs of the folders containing
# fastq files. To add a new study to this process, a row needs to be added for
# that study in this CSV file and any additional fixes to that study's metadata
# should be added in the script below.

library(synapser)
library(stringr)
library(dplyr)

# Synapse IDs used in this script that are not in Model_AD_SynID_list.csv
syn_metadata_folder_id <- "syn61850200"

synLogin(silent = TRUE)
tmp_dir <- file.path("data", "tmp")
dir.create(tmp_dir, showWarnings = FALSE)

metadata_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"))

# Process each study's metadata ------------------------------------------------

metadata <- lapply(1:nrow(metadata_list), function(N) {
  row <- metadata_list[N, ]
  print(row$Study)

  ## Fetch metadata files ------------------------------------------------------

  assay_file <- synGet(row$Metadata_Assay,
                       downloadLocation = tmp_dir,
                       ifcollision = "overwrite.local")
  biospec_file <- synGet(row$Metadata_Biospecimen,
                         downloadLocation = tmp_dir,
                         ifcollision = "overwrite.local")
  individual_file <- synGet(row$Metadata_Individual,
                            downloadLocation = tmp_dir,
                            ifcollision = "overwrite.local")

  assay_df <- read.csv(assay_file$path)
  biospec_df <- read.csv(biospec_file$path)
  individual_df <- read.csv(individual_file$path)


  ## Study-specific fixes ------------------------------------------------------

  # Jax.IU.Pitt_APOE4.Trem2.R47H: We are temporarily using a file from the
  # staging version of the project, not what is released in the portal. This
  # file isn't technically in the format of a portal assay file but has the
  # information we need, with some extra restructuring needed. This file has two
  # rows for some specimenIDs because of sequencing on two separate lanes. We
  # remove duplicate rows since the fastqs we are using have both lanes merged.
  if (row$Study == "Jax.IU.Pitt_APOE4.Trem2.R47H") {
    assay_df <- assay_df %>%
      dplyr::rename(individualID = animalName,
                    specimenID = sampleName) %>%

      # Make totalReads a number
      mutate(totalReads = str_replace_all(totalReads, ",", ""),
             totalReads = as.numeric(totalReads)) %>%

      # Remove references to specific fastq files
      select(!contains("fastq"), -path, -study, -lane) %>%

      # Add totalReads for both lanes together. Some specimenIDs were sequenced
      # in two different batches so we combine the batch names
      group_by_at(setdiff(colnames(.), c("totalReads", "sequencingBatch"))) %>%
      summarize(totalReads = sum(totalReads),
                sequencingBatch = paste(sort(unique(sequencingBatch)), collapse = ";"),
                .groups = "drop") %>%

      # Add missing columns to make the format match assay metadata files
      mutate(platform = NA, RIN = NA,
             rnaBatch = sequencingBatch,
             libraryBatch = sequencingBatch)

  } else if (row$Study == "UCI_5XFAD") {
    # UCI_5XFAD: The specimen IDs in the assay metadata do not match what is in
    # the biospecimen metadata, so we need to fix them. In the assay metadata
    # they are of the format "<individualID>(H or C)_RNAseq", (e.g.
    # "295C_RNAseq") while in the biospecimen metadata they are of the format
    # "<individualID>r(h or c)", (e.g. "295rc").
    assay_df$specimenID <- str_replace(assay_df$specimenID, "_RNAseq", "")
    assay_df$specimenID <- str_replace(assay_df$specimenID, "H", "rh")
    assay_df$specimenID <- str_replace(assay_df$specimenID, "C", "rc")

  } else if (row$Study == "UCI_hAbeta_KI") {
    # UCI_hAbeta_KI: Some platform entries are mis-labeled as NextSeq501,
    # NextSeq502, or NextSeq503 when they should all be NextSeq500. The
    # biospecimen metadata file is missing the "samplingAge" column.
    assay_df$platform <- str_replace(assay_df$platform,
                                     "NextSeq50[1|2|3]",
                                     "NextSeq500")
    biospec_df$samplingAge <- NA
  }

  # Studies Jax.IU.Pitt_5XFAD, UCI_3xTg-AD, UCI_ABCA7, UCI_PrimaryScreen,
  # UCI_Trem2_Cuprizone, and UCI_Trem2-R47H_NSS need no specialized corrections
  # in this section


  ## Merge all metadata together -----------------------------------------------

  # Only keep rows that exist in all dataframes so we only retain RNA
  # seq-related samples.
  combined_df <- merge(assay_df, biospec_df, all = FALSE) %>%
    merge(individual_df, all = FALSE)

  combined_df$study_name <- row$Study


  ## Fix some specimenIDs post-merge -------------------------------------------

  # UCI_hAbeta_KI has parenthesis in some specimenIDs, which we remove to avoid
  # issues downstream
  combined_df$specimenID <- str_replace(combined_df$specimenID, "\\(", "")
  combined_df$specimenID <- str_replace(combined_df$specimenID, "\\)", "")

  # UCI_hAbeta_KI has commas in some individualIDs and specimenIDs. The commas
  # will cause problems with CSV files downstream so we remove them here.
  combined_df$specimenID <- str_replace(combined_df$specimenID, "\\,", "")
  combined_df$individualID <- str_replace(combined_df$individualID, "\\,", "")

  # What the specimen IDs will look like when the counts files are read in by
  # R, to make it easier to match to the metadata. "R_safe_specimenID" is
  # for reading in individual study files, while "merged_file_specimenID" is
  # what they will look like in the merged counts files containing all studies.
  combined_df$R_safe_specimenID <- make.names(combined_df$specimenID)
  combined_df$merged_file_specimenID <- make.names(
    paste(combined_df$study_name, combined_df$specimenID, sep = ".")
  )

  # General fix to all studies, we need the "treatmentType" and "treatmentDose"
  # columns for the Cuprizone study but some of the other studies don't have
  # these fields
  if (!hasName(combined_df, "treatmentType")) {
    combined_df$treatmentType <- NA
  }
  if (!hasName(combined_df, "treatmentDose")) {
    combined_df$treatmentDose <- NA
  }


  ## Create final study-specific data frame ------------------------------------

  # Filter to columns of interest and ensure that all columns are in the same
  # order across all studies. Also make sure all rows are unique (some studies
  # have duplicate rows in one or more files).
  combined_df <- combined_df %>%
    select(individualID, specimenID, R_safe_specimenID, merged_file_specimenID,
           platform, RIN, rnaBatch, libraryBatch, sequencingBatch, organ,
           tissue, samplingAge, treatmentType, treatmentDose, sex, ageDeath,
           ageDeathUnits, genotype, genotypeBackground, modelSystemName,
           study_name) %>%
    distinct()
  return(combined_df)
})

# Bind into one data frame
metadata_combined <- do.call(rbind, metadata)


# Fixes to genotype names ------------------------------------------------------

# Standardize all genotypes to the MODEL-AD approved values.
# Note: "C57BL6J" doesn't exist on the approved list but I think it's ok to
#       leave it as-is since it's used as a general control for multiple
#       genotypes.
# Note: "Trem2-KO" and "Trem2-R47H_CSS_homozygous" don't exist on the approved
#       list but should be on it, so we leave these values as-is.
geno_map <- c("3xTg-AD_homozygous" = "3xTg-AD_carrier",
              "3XTg-AD_noncarrier" = "3xTg-AD_noncarrier",
              "5XFAD_hemizygous" = "5XFAD_carrier",
              "Abca7_V1599M_homozygous" = "Abca7-V1599M_homozygous",
              "Abca7V1599M_noncarrier" = "Abca7-V1599M_WT",
              "ABI3_S209F_homozygous" = "Abi3-S209F_homozygous",
              "ABI3_S209F_noncarrier" = "Abi3-S209F_WT",
              "BIN1_K358R_noncarrier" = "Bin1-K358R_WT",
              "hABKI  HO" = "hAbeta-KI_LoxP_homozygous",
              "hABKI  WT" = "hAbeta-KI_LoxP_WT",
              "Homozygous" = "homozygous",
              "Heterozygous" = "heterozygous",
              "NA; NA" = NA,
              "PICALM_H458R_homozygous" = "Picalm-H458R_homozygous",
              "PICALM_H458R_noncarrier" = "Picalm-H458R_WT",
              "SPI1_homozygous" = "Spi1-rs1377416_homozygous",
              "SPI1_noncarrier" = "Spi1-rs1377416_WT",
              "Trem2_R47H_homozygous" = "Trem2-R47H_homozygous",
              "Trem2_R47H_noncarrier" = "Trem2-R47H_WT",
              "TREM2R47H_heterozygous" = "Trem2-R47H_heterozygous",
              "TREM2R47H_homozygous" = "Trem2-R47H_homozygous",
              "TREM2R47H_noncarrier" = "Trem2-R47H_WT"
              )

for (G in 1:length(geno_map)) {
  metadata_combined$genotype <- str_replace_all(metadata_combined$genotype,
                                                names(geno_map)[G],
                                                geno_map[G])
}

# Genotypes should be semicolon-separated, not comma-separated
metadata_combined$genotype <- str_replace(metadata_combined$genotype, ",", ";")


# Save to file and upload to Synapse -------------------------------------------

write.csv(metadata_combined, file.path("data", "Model_AD_merged_metadata.csv"),
          row.names = FALSE, quote = FALSE)

syn_file <- File(file.path("data", "Model_AD_merged_metadata.csv"),
                 parent = syn_metadata_folder_id)

all_syn_ids <- c(metadata_list$Metadata_Assay,
                 metadata_list$Metadata_Biospecimen,
                 metadata_list$Metadata_Individual)
github_link <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/01_Harmonize_Metadata.R"

synStore(syn_file,
         used = all_syn_ids,
         executed = github_link,
         forceVersion = FALSE)
