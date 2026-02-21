# This script pulls in all metadata files necessary for RNA seq analysis from
# all studies included in the analysis and cleans up a few fields. Generally 3
# metadata files are needed: an individual metadata file, a biospecimen metadata
# file, and an RNA seq assay metadata file, which are merged together so that
# the final data frame only contains rows for specimens that have RNA seq data.
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
#   UCI_Trem2-R47H_NSS
#
# The CSV file "Model_AD_SynID_list.csv" was created by hand and lists the
# Synapse IDs of all metadata files, some of which do not currently exist on the
# portal and cannot be queried for. To add a new study to this process, a row
# needs to be added for that study in this CSV file and any additional fixes to
# that study's metadata should be added in the script below.

library(synapser)
library(stringr)
library(dplyr)

folder_syn_ids <- config::get("folder_syn_ids", config = "default")
studies <- config::get("studies", config = "default")

synLogin(silent = TRUE)
tmp_dir <- file.path("output", "tmp")
dir.create(tmp_dir, showWarnings = FALSE)

metadata_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"),
                          comment.char = "#") |>
  subset(Study %in% studies)

stopifnot(all(studies %in% metadata_list$Study))


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

  for (syn_file in c(assay_file, biospec_file, individual_file)) {
    if (!syn_file$isLatestVersion) {
      warning(
        str_glue("'{syn_file$name}' ({syn_file$id} v{syn_file$versionNumber}) ",
                 "is not the most recent version on Synapse! \nThe specified ",
                 "version ({syn_file$versionNumber}) will be used for ",
                 "harmonization unless 'Model_AD_SynIDlist.csv' is updated ",
                 "with with a newer version.")
      )
    }
  }

  assay_df <- read.csv(assay_file$path)
  biospec_df <- read.csv(biospec_file$path)
  individual_df <- read.csv(individual_file$path)


  ## Study-specific fixes ------------------------------------------------------

  if (row$Study == "Jax.IU.Pitt_LOAD2.PrimaryScreen") {
    # Genotypes in this study leave off the "LOAD2" at the front and it needs to
    # be added. We make this modification here instead of below with the rest of
    # the genotype mapping so we don't affect genotypes from other studies.
    geno_mods <- which(!grepl("LOAD2", individual_df$genotype) &
                         !(individual_df$genotype %in% c("WT", "", "NA_Inconclusive")))
    individual_df$genotype[geno_mods] <- paste0("LOAD2.",
                                                individual_df$genotype[geno_mods])

    # Change "WT" to "C57BL6J"
    individual_df$genotype[individual_df$genotype == "WT"] <- "C57BL6J"

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

  # Studies Jax.IU.Pitt_5XFAD, Jax.IU.Pitt_APOE4.Trem2.R47H, UCI_3xTg-AD,
  # UCI_ABCA7, UCI_Bin1K358R, UCI_Clu-h2kbKI, UCI_PrimaryScreen,
  # UCI_Trem2_Cuprizone, and UCI_Trem2-R47H_NSS need no specialized corrections
  # in this section

  # General corrections -- remove any rows with NA or "" IDs. Some files have
  # empty rows at the end.
  assay_df <- subset(assay_df,
                     !is.na(specimenID) & nchar(specimenID) > 0)
  biospec_df <- subset(biospec_df,
                       !is.na(specimenID) & nchar(specimenID) > 0 &
                         !is.na(individualID) & nchar(individualID) > 0)
  individual_df <- subset(individual_df,
                          !is.na(individualID) & nchar(individualID) > 0)

  # Check for issues with specimenID / individualID in each file.
  # If this section produces errors, the study needs corrections added above.
  stopifnot(all(!is.na(assay_df$specimenID)) & all(nchar(assay_df$specimenID) > 0))
  stopifnot(!any(duplicated(assay_df$specimenID)))
  stopifnot(all(assay_df$specimenID %in% biospec_df$specimenID))

  stopifnot(all(!is.na(biospec_df$specimenID)) & all(nchar(biospec_df$specimenID) > 0))
  stopifnot(all(!is.na(biospec_df$individualID)) & all(nchar(biospec_df$individualID) > 0))
  stopifnot(all(biospec_df$individualID %in% individual_df$individualID))

  stopifnot(all(!is.na(individual_df$individualID)) & all(nchar(individual_df$individualID) > 0))
  stopifnot(!any(duplicated(individual_df$individualID)))


  ## Merge all metadata together -----------------------------------------------

  # Only keep rows that exist in all data frames so we only retain RNA
  # seq-related samples.
  combined_df <- merge(assay_df, biospec_df, all = FALSE) |>
    merge(individual_df, all = FALSE)

  combined_df$study <- row$Study


  ## Fix some specimenIDs post-merge -------------------------------------------

  # What the specimen IDs will look like when the counts files are read in by
  # R, to make it easier to match to the metadata. "R_safe_specimenID" is
  # for reading in individual study files, while "unique_specimenID" is
  # a specimenID guaranteed to be unique across all studies, since some studies
  # share some specimenIDs.
  combined_df$R_safe_specimenID <- make.names(combined_df$specimenID)
  combined_df$unique_specimenID <- make.names(
    paste(combined_df$study, combined_df$specimenID, sep = ".")
  )


  ## Create final study-specific data frame ------------------------------------

  # Filter to columns of interest and ensure that all columns are in the same
  # order across all studies. Also make sure all rows are unique (some studies
  # have duplicate rows in one or more files).
  combined_df <- combined_df |>
    select(
      # individual metadata
      study, individualID, specimenID, R_safe_specimenID, unique_specimenID,
      sex, ageDeath, ageDeathUnits, genotype, genotypeBackground,
      # biospec metadata
      tissue,
      # assay metadata
      platform, RIN, rnaBatch, libraryBatch, sequencingBatch
    ) |>
    distinct()
  return(combined_df)
})

# Bind into one data frame so we can make genotype changes all at once
metadata_combined <- do.call(rbind, metadata)


# Fixes to genotype names ------------------------------------------------------

# Standardize all genotypes to the MODEL-AD approved values.
geno_map <- c(
  "3xTg-AD_homozygous" = "3xTg-AD_carrier",
  "3XTg-AD_noncarrier" = "3xTg-AD_noncarrier",
  "5XFAD_hemizygous" = "5XFAD_carrier",
  "Homozygous" = "homozygous",
  "Heterozygous" = "heterozygous",
  "NA_Inconclusive" = NA,
  # Not relevant to set of studies used for analysis but saving for
  # reproducibility.
  "Abca7_V1599M_homozygous" = "Abca7-V1599M_homozygous",
  "Abca7V1599M_noncarrier" = "Abca7-V1599M_WT",
  "ABI3_S209F_homozygous" = "Abi3-S209F_homozygous",
  "ABI3_S209F_noncarrier" = "Abi3-S209F_WT",
  "BIN1_K358R_noncarrier" = "Bin1-K358R_WT",
  "hABKI  HO" = "hAbeta-KI_LoxP_homozygous",
  "hABKI  WT" = "hAbeta-KI_LoxP_WT",
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

# Standardize genotypeBackground values
metadata_combined <- metadata_combined |>
  mutate(genotypeBackground = case_match(
    genotypeBackground,
    c("B6", "C57BL6J")  ~ "C57BL/6J",
    .default = genotypeBackground
  ))


# TODO the genotype names for LOAD2 Primary Screen may not be in the right format


# Save to files and upload to Synapse ------------------------------------------

dir.create(file.path("output", "metadata"), showWarnings = FALSE)

# Split the metadata up by study and write to individual files
for (study_name in unique(metadata_combined$study)) {
  study_data <- subset(metadata_combined, study == study_name)
  study_file <- file.path("output", "metadata",
                          paste0(study_name, "_harmonized_metadata.csv"))

  write.csv(study_data, study_file, row.names = FALSE, quote = FALSE)

  syn_file <- File(study_file, parent = folder_syn_ids$metadata)

  all_syn_ids <- subset(metadata_list, Study == study_name) |>
    select(-Study) |>
    as.character()

  github_link <- paste0(config::get("github_repo_url", config = "default"),
                        "/blob/main/02_Harmonize_Metadata.R")

  synStore(syn_file,
           used = all_syn_ids,
           executed = github_link,
           forceVersion = FALSE)
}
