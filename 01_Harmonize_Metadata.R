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

synLogin()
tmp_dir <- file.path("data", "tmp")
dir.create(tmp_dir, showWarnings = FALSE)

metadata_list <- read.csv(file.path("data", "Model_AD_SynID_list.csv"))

# Process each study's metadata ------------------------------------------------

metadata <- lapply(1:nrow(metadata_list), function(N) {
  row <- metadata_list[N, ]
  print(row$Study)

  ## Fetch metadata files ------------------------------------------------------

  assay_file <- synGet(row$Metadata_Assay,
                       downloadLocation = tmp_dir)
  biospec_file <- synGet(row$Metadata_Biospecimen,
                         downloadLocation = tmp_dir)
  individual_file <- synGet(row$Metadata_Individual,
                            downloadLocation = tmp_dir)

  assay_df <- read.csv(assay_file$path)
  biospec_df <- read.csv(biospec_file$path)
  individual_df <- read.csv(individual_file$path)


  ## Study-specific fixes ------------------------------------------------------

  # Jax.IU.Pitt_5XFAD: the "ageDeath" field is missing in the individual
  # metadata, so we calculate it from dateDeath and dateBirth. There are some
  # NA values in these fields but not in any individuals with RNA seq data.
  if (row$Study == "Jax.IU.Pitt_5XFAD") {
    ageDeath <- as.Date(individual_df$dateDeath, format = "%m/%d/%Y") -
      as.Date(individual_df$dateBirth, format = "%m/%d/%Y")
    individual_df$ageDeath <- as.numeric(ageDeath) / 12 # days -> months
    individual_df$ageDeathUnits <- "months"
  }

  # Jax.IU.Pitt_APOE4.Trem2.R47H: For some reason this study doesn't have any of
  # the brain samples in its assay file so we can't merge those files. Instead
  # we create a fake dataframe with the expected structure of an assay metadata
  # file.
  if (row$Study == "Jax.IU.Pitt_APOE4.Trem2.R47H") {
    # The %in% is necessary instead of != in this statement due to NAs
    brain_only <- subset(biospec_df, organ == "brain" &
                           !(specimenIdSource %in% c("soluble")))

    # Create a fake assay metadata dataframe with empty values for fields that
    # need to be there but we don't have information for.
    assay_df <- data.frame(specimenID = brain_only$specimenID,
                           platform = NA,
                           RIN = NA,
                           rnaBatch = NA,
                           libraryBatch = NA,
                           sequencingBatch = NA)
  }

  # UCI_3xTg-AD and UCI_5XFAD: some capitalization differences in the individual
  # metadata file
  if (row$Study %in% c("UCI_3xTg-AD", "UCI_5XFAD")) {
    individual_df$stockNumber <- individual_df$StockNumber
    individual_df$officialName <- individual_df$OfficialName
  }

  # UCI_5XFAD: The specimen IDs in the assay metadata do not match what is in
  # the biospecimen metadata, so we need to fix them. In the assay metadata
  # they are of the format "<individualID>(H or C)_RNAseq", (e.g. "295C_RNAseq")
  # while in the biospecimen metadata they are of the format
  # "<individualID>r(h or c)", (e.g. "295rc").
  #
  if (row$Study == "UCI_5XFAD") {
    assay_df$specimenID <- str_replace(assay_df$specimenID, "_RNAseq", "")
    assay_df$specimenID <- str_replace(assay_df$specimenID, "H", "rh")
    assay_df$specimenID <- str_replace(assay_df$specimenID, "C", "rc")
  }

  # UCI_hAbeta_KI: Some platform entries are mis-labeled as NextSeq501,
  # NextSeq502, or NextSeq503 when they should all be NextSeq500. The
  # biospecimen metadata file is missing the "samplingAge" column.
  if (row$Study == "UCI_hAbeta_KI") {
    assay_df$platform <- str_replace(assay_df$platform,
                                     "NextSeq50[1|2|3]",
                                     "NextSeq500")
    biospec_df$samplingAge <- NA
  }

  # Studies UCI_ABCA7, UCI_PrimaryScreen, UCI_Trem2_Cuprizone, and
  # UCI_Trem2-R47H_NSS need no corrections in this section


  ## Merge all metadata together -----------------------------------------------

  # Only keep rows that exist inall dataframes so we only retain RNA seq-related
  # samples.
  combined_df <- merge(assay_df, biospec_df,
                       by = "specimenID",
                       all = FALSE)

  combined_df <- merge(combined_df, individual_df,
                       by = "individualID",
                       all = FALSE)

  combined_df$study_name <- row$Study

  # Filter to columns of interest and ensure that all columns are in the same
  # order across all studies. Also make sure all rows are unique (some studies
  # have duplicate rows in one or more files).
  combined_df <- select(combined_df, individualID, specimenID, platform, RIN,
                        rnaBatch, libraryBatch, sequencingBatch, organ, tissue,
                        samplingAge, sex, ageDeath, ageDeathUnits, genotype,
                        genotypeBackground, stockNumber, officialName,
                        modelSystemName, study_name) %>%
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
# Note: "APOE4-KI_noncarrier" is listed as "APOE4-KI_WT" in the approved list,
#       however this is a transgene so it seems more appropriate to label it
#       as "noncarrier", so I've done that here.
geno_map <- c("3xTg-AD_homozygous" = "3xTg-AD_carrier",
              "3XTg-AD_noncarrier" = "3xTg-AD_noncarrier",
              "5XFAD_hemizygous" = "5XFAD_carrier",
              "Abca7_V1599M_homozygous" = "Abca7-V1599M_homozygous",
              "Abca7V1599M_noncarrier" = "Abca7-V1599M_WT",
              "ABI3_S209F_homozygous" = "Abi3-S209F_homozygous",
              "ABI3_S209F_noncarrier" = "Abi3-S209F_WT",
              "APOE4_heterozygous" = "APOE4-KI_heterozygous",
              "APOE4_homozygous" = "APOE4-KI_homozygous",
              "APOE4_noncarrier" = "APOE4-KI_noncarrier",
              "BIN1_K358R_noncarrier" = "Bin1-K358R_WT",
              "hABKI  HO" = "hAbeta-KI_LoxP_homozygous",
              "hABKI  WT" = "hAbeta-KI_LoxP_WT",
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
  metadata_combined$genotype <- str_replace(metadata_combined$genotype,
                                            names(geno_map)[G],
                                            geno_map[G])
}

# Genotypes should be semicolon-separated, not comma-separated
metadata_combined$genotype <- str_replace(metadata_combined$genotype, ",", ";")


# Study-specific fixes post-merge ----------------------------------------------

# UCI_hAbeta_KI: has parenthesis in some specimenIDs, which we remove to avoid
# issues downstream
metadata_combined$specimenID <- str_replace(metadata_combined$specimenID, "\\(", "")
metadata_combined$specimenID <- str_replace(metadata_combined$specimenID, "\\)", "")

# Jax.IU.Pitt_5XFAD and UCI_hAbeta_KI: have commas in some individualIDs and
# specimenIDs. The commas will cause problems with CSV files downstream so we
# remove them here.
metadata_combined$specimenID <- str_replace(metadata_combined$specimenID, "\\,", "")
metadata_combined$individualID <- str_replace(metadata_combined$individualID, "\\,", "")


# Save to file and upload to Synapse -------------------------------------------

write.csv(metadata_combined, file.path("data", "Model_AD_merged_metadata.csv"),
          row.names = FALSE, quote = FALSE)

syn_file <- File(file.path("data", "Model_AD_merged_metadata.csv"),
                 parent = "syn61850200")

all_syn_ids <- c(metadata_list$Metadata_Assay,
                 metadata_list$Metadata_Biospecimen,
                 metadata_list$Metadata_Individual)
github_link <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/Harmonize_Metadata.R"

synStore(syn_file,
         used = all_syn_ids,
         executed = github_link,
         forceVersion = FALSE)
