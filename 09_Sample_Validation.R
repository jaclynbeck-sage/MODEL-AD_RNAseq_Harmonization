library(synapser)
library(synapserutils)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

source("validation_functions.R")

# Set up -----------------------------------------------------------------------

file_syn_ids <- config::get("file_syn_ids", config = "default")
folder_syn_ids <- config::get("folder_syn_ids", config = "default")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url", config = "default"),
                 "/blob/main/09_Sample_Validation.Rmd")
tmp_dir <- file.path("output", "tmp")


# Utility functions ------------------------------------------------------------

# Get a data frame with specific genes, merged with metadata
make_counts_df <- function(metadata, counts, symbol_map, genes) {
  genes_sub <- subset(symbol_map, gene_symbol %in% genes)
  rownames(genes_sub) <- genes_sub$ensembl_gene_id

  counts_sub <- counts[genes_sub$ensembl_gene_id, ]
  rownames(counts_sub) <- genes_sub[rownames(counts_sub), "gene_symbol"]

  counts_df <- as.data.frame(t(counts_sub)) %>%
    merge(metadata, by.x = "row.names", by.y = "merged_file_specimenID") %>%
    dplyr::rename(merged_file_specimenID = Row.names)
}


score_quality <- function(quality) {
  ifelse(is.na(quality), "deletion",
         ifelse(quality < 20, "low",
                ifelse(quality < 100, "moderate",
                       "high")))
}


get_variant_mismatches <- function(metadata, geno_info, genotype_pattern,
                                   mutation_match = NULL,
                                   gene_symbol_match = NULL,
                                   total_positions = 0,
                                   rm_na = TRUE) {
  carrier_samples <- subset(metadata, grepl(genotype_pattern, genotype))

  # Get every sample in the same studies as the carrier samples
  study_samples <- subset(metadata, study_name %in% carrier_samples$study_name)

  geno_sub <- subset(geno_info, merged_file_specimenID %in% study_samples$merged_file_specimenID &
                       (gene_symbol %in% gene_symbol_match |
                          mutation %in% mutation_match)) %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype, gene_symbol) %>%
    summarize(count = n(),
              avg_quality = mean(quality, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(evidence = score_quality(avg_quality),
           is_carrier = grepl(genotype_pattern, genotype))

  # Add in samples that don't show up in geno_sub
  missing_samples <- study_samples %>%
    subset(!(merged_file_specimenID %in% geno_sub$merged_file_specimenID)) %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype) %>%
    reframe(gene_symbol = unique(geno_sub$gene_symbol),
            count = 0, avg_quality = 0, evidence = "none",
            is_carrier = grepl(genotype_pattern, genotype))

  geno_sub <- rbind(geno_sub, missing_samples)

  geno_final <- geno_sub %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype, is_carrier) %>%
    summarize(total_found = sum(count),
              evidence = paste(unique(evidence), collapse = ","),
              .groups = "drop")

  if (rm_na) {
    geno_final$total_found[geno_final$evidence == "deletion"] <- 0
  }

  geno_final <- geno_final %>%
    mutate(est_genotype =
             ifelse(total_found >= total_positions, "carrier",
                    ifelse(total_found == 0 & evidence == "none", "noncarrier",
                           "ambiguous")))

  return(geno_final)
}


print_variant_mismatches <- function(var_mismatches, genotype_name) {
  wt_mismatches <- subset(var_mismatches, !is_carrier)
  carrier_mismatches <- subset(var_mismatches, is_carrier)

  if (nrow(wt_mismatches) > 0) {
    print(paste("Detected", genotype_name, "variants in",
                nrow(wt_mismatches), "WT samples:"))

    print(select(wt_mismatches, study_name, specimenID, genotype, est_genotype))
  } else {
    print(paste("No WT samples have detected", genotype_name, "variants."))
  }

  print("")

  if (nrow(carrier_mismatches) > 0) {
    print(paste(nrow(carrier_mismatches), "carrier samples",
                "had no detected", genotype_name, "variants:"))
    print(select(carrier_mismatches, study_name, specimenID, genotype, est_genotype))
  } else {
    print(paste("No carrier samples are missing detected", genotype_name,
                "variants."))
  }
}


print_expression_mismatches <- function(mismatch_df, genotype_name) {
  if (nrow(mismatch_df) > 0) {
    print(paste(nrow(mismatch_df), "samples have a", genotype_name,
                "expression mismatch:"))
    print(select(mismatch_df, -merged_file_specimenID))
  } else {
    print(paste("No samples have a", genotype_name, "expression mismatch."))
  }
}


# Load counts and metadata -----------------------------------------------------

metadata_file <- synGet(file_syn_ids$merged_metadata,
                        downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")
counts_file <- synGet(file_syn_ids$merged_counts,
                      downloadLocation = tmp_dir,
                      ifcollision = "overwrite.local")
symbol_map_file <- synGet(file_syn_ids$symbol_map,
                          downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

metadata_all <- read.csv(metadata_file$path)
counts <- read.delim(counts_file$path, header = TRUE, row.names = 1) %>%
  # Get rid of transcript_id column
  select(-transcript_id.s.)

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, merged_file_specimenID %in% colnames(counts)) %>%
  select(individualID, specimenID, merged_file_specimenID, sex, genotype, ageDeath, study_name)
counts <- counts[, metadata_all$merged_file_specimenID]

symbol_map <- read.csv(symbol_map_file$path)

# Convert counts to CPM so samples are comparable to each other
lib_sizes <- colSums(counts)
counts <- sweep(counts, 2, lib_sizes, "/") * 1e6


# Find vcf files on Synapse ----------------------------------------------------

# TODO: ideally these would be annotated with the specimenID / study name but
# for now we walk the folder structure and use the folder name for each sample
# to determine the specimenID

studies <- synGetChildren(folder_syn_ids$genotype_validation,
                          includeTypes = list("folder"))$asList()
names(studies) <- sapply(studies, "[[", "name")

studies <- studies[names(studies) %in% metadata_all$study_name]

study_vcfs <- lapply(studies, function(study) {
  folder_structure <- synapserutils::walk(study$id, includeTypes = list("file"))$asList()

  samples_df <- lapply(folder_structure, function(item) {
    # The structure of each item that contains an actual file is a list:
    #   item[[1]] = list(1 = path of containing folder, 2 = synapse id of folder)
    #   item[[2]] = list(empty), or a list of sub-folders
    #   item[[3]] = list of files, where each item is list(1 = file name, 2 = synapse id)
    if (length(item) < 3 || lengths(item)[3] == 0) {
      return(NULL)
    }

    folder_name <- pluck(item, 1, 1)
    specimenID <- basename(folder_name)

    # Find the file that ends in .gz, not .gz.tbi
    vcf_filenames <- sapply(item[[3]], "[[", 1)
    ind <- which(grepl("gz$", vcf_filenames))
    return(data.frame(study_name = study$name,
                      specimenID = specimenID,
                      filename = vcf_filenames[ind],
                      syn_id = pluck(item, 3, ind, 2)))
  })

  samples_df <- do.call(rbind, samples_df)
  return(samples_df)
})


# Download and parse vcf files -------------------------------------------------

# Read the intervals file used to call variants and turn it into one row per
# genome position rather than one row per range of positions.
# The intervals file coordinates are -1 and +1 off from the real coordinate(s)
# as reported in the vcf file:
#   start = real coordinate - 1 and
#   end = real coordinate + 1
# We undo this operation so the values in the "position" column match the values
# in the vcf files.
intervals_file <- synGet(file_syn_ids$intervals,
                         downloadLocation = tmp_dir,
                         ifcollision = "overwrite.local")

intervals <- read.delim(intervals_file$path, header = FALSE) %>%
  dplyr::rename(chromosome = V1, start = V2, end = V3, mutation = V4) %>%
  rowwise() %>%
  # adds as many new rows as necessary to cover the full range of positions
  reframe(chromosome = chromosome,
          mutation = mutation,
          position = (start+1):(end-1))

geno_info <- lapply(study_vcfs, function(study_df) {
  mutation_df <- apply(study_df, 1, function(row) {
    vcf_file <- synGet(row[["syn_id"]], downloadLocation = tmp_dir,
                       ifcollision = "overwrite.local")

    vcf <- NULL

    # If the sample didn't have any detected mutations, the vcf file will
    # contain only comments and no data, which will cause the read.delim
    # function to throw an error.
    try({
      vcf <- read.delim(gzfile(vcf_file$path), header = FALSE, comment.char = "#") %>%
        dplyr::rename(chromosome = V1, position = V2, ref = V4, alt = V5,
                      quality = V6) %>%
        select(chromosome, position, ref, alt, quality) %>%
        mutate(study_name = row[["study_name"]],
               specimenID = row[["specimenID"]]) %>%
        merge(intervals)
    }, silent = TRUE)

    # If no mutations found, this will be NULL
    return(vcf)
  })

  mutation_df <- do.call(rbind, mutation_df) %>%
    mutate(gene_symbol = str_replace(mutation, "[-\\*].*", ""))

  return(mutation_df)
})

geno_info <- do.call(rbind, geno_info)
geno_info$quality[geno_info$quality == "."] <- NA
geno_info$quality <- as.numeric(geno_info$quality)
geno_info <- merge(geno_info, metadata_all) %>%
  mutate(evidence = score_quality(quality))


# Set up the valid samples data frame ------------------------------------------

#valid_samples <- metadata_all %>%
#  select(individualID, specimenID, merged_file_specimenID, study_name)
valid_samples <- data.frame()


# Validation of mouse sex ------------------------------------------------------

xy_df <- make_counts_df(metadata_all, counts, symbol_map,
                        c("Xist", "Eif2s3y", "Ddx3y")) %>%
  # More than 10 CPM for Xist -- there is low-level expression even in males
  mutate(est_female = Xist > 10,
         # Expression of Y-related genes should be zero for females, so we
         # use anything > 1 CPM for males and assume < 1 CPM might be noise
         est_male = Eif2s3y > 1 | Ddx3y > 1)

# Mark valid/invalid conditions
xy_df$valid_sex <- FALSE
xy_df$valid_sex[xy_df$sex == "female" &
                  xy_df$est_female == TRUE &
                  xy_df$est_male == FALSE] <- TRUE
xy_df$valid_sex[xy_df$sex == "male" &
                  xy_df$est_female == FALSE &
                  xy_df$est_male == TRUE] <- TRUE

# Mark valid samples in the data frame
sex_matches <- select(xy_df, merged_file_specimenID, valid_sex)
#valid_samples <- merge(valid_samples, sex_matches, all = TRUE)

# For printing
mismatches <- subset(xy_df, !valid_sex) %>%
  select(study_name, individualID, specimenID, sex, est_female, est_male, Xist, Eif2s3y, Ddx3y) %>%
  dplyr::rename(reported_sex = sex)
mismatches$est_sex <- "unknown"
mismatches$est_sex[mismatches$est_female & !mismatches$est_male] <- "female"
mismatches$est_sex[mismatches$est_male & !mismatches$est_female] <- "male"

print(paste(nrow(mismatches), "samples have mismatched sex:"))
print(select(mismatches, study_name, specimenID, reported_sex, est_sex,
             Xist, Eif2s3y, Ddx3y))

# TODO 4 LOAD2 samples fail but gene expression indicates they are correct. It just
# falls outside of the thresholds I set

# Validation of mouse genotype -------------------------------------------------

## Jax.IU.Pitt_5XFAD and UCI_5XFAD ---------------------------------------------

meta_5x <- subset(metadata_all,
                  study_name %in% c("Jax.IU.Pitt_5XFAD", "UCI_5XFAD"))

valid_5x <- validate_5x(meta_5x, geno_info, counts, symbol_map)

valid_samples <- rbind(valid_samples, valid_5x$valid)

# Tmp PCA

counts_5x <- log2(counts[, meta_5x$merged_file_specimenID] + 0.5)

for (study in unique(meta_5x$study_name)) {
  meta_study <- subset(meta_5x, study_name == study) |>
    merge(valid_5x$detail)
  vars <- matrixStats::rowVars(as.matrix(counts_5x[, meta_study$merged_file_specimenID]))
  hv <- names(sort(vars, decreasing = TRUE))[1:2000]
  pc <- prcomp(t(counts_5x[hv, meta_study$merged_file_specimenID]))
  pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                   by.y = "merged_file_specimenID") %>%
    dplyr::rename(merged_file_specimenID = Row.names) %>%
    mutate(
      has_5x = grepl("5XFAD_carrier", genotype),
      mismatch_type = case_when(
        valid_5x_variant & valid_5x_expression ~ "none",
        valid_5x_variant & !valid_5x_expression ~ "expression",
        !valid_5x_variant & valid_5x_expression ~ "variant",
        !valid_5x_variant & !valid_5x_expression ~ "expression + variant"
      )
    )
  plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = has_5x)) +
    geom_point() + ggtitle(study)
  print(plt)
}

# TODO GT19_12887 (non-carrier) from Jax 5X study has all 6 variants detected
# but expresses only a small amount of APP (11.28 CPM) and PSEN1 (0.85 CPM)
# and clusters with other non-carriers in a PCA


## Jax.IU.Pitt_APOE4.Trem2.R47H ------------------------------------------------

# APOE4 KI can be validated with gene expression only, and Trem2-R47H can be
# validated with variant calling.
meta_load1 <- subset(metadata_all, study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H")

valid_apoe <- validate_APOE4_KI(meta_load1, counts, symbol_map)
valid_trem2 <- validate_Trem2_R47H(meta_load1, geno_info)

valid_all <- merge(valid_apoe$valid, valid_trem2$valid,
                   by = "merged_file_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(merged_file_specimenID, valid)

valid_samples <- rbind(valid_samples, valid_all)

# Tmp PCA

counts_apoe <- log2(counts[, meta_load1$merged_file_specimenID] + 0.5)

meta_study <- subset(meta_load1, !is.na(genotype)) |>
  merge(valid_apoe$detail) |>
  merge(valid_trem2$detail)

vars <- matrixStats::rowVars(as.matrix(counts_apoe[, meta_study$merged_file_specimenID]))
hv <- names(sort(vars, decreasing = TRUE))[1:2000]

pc <- prcomp(t(counts_apoe[hv, meta_study$merged_file_specimenID]))

pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                 by.y = "merged_file_specimenID") |>
  dplyr::rename(merged_file_specimenID = Row.names) |>
  mutate(
    is_load1 = grepl("APOE4-KI_(homo|hetero)", genotype),
    mismatch_type = case_when(
      valid_apoe4_expression & valid_trem2_r47h_variant ~ "none",
      !valid_apoe4_expression & valid_trem2_r47h_variant ~ "APOE4 expression",
      valid_apoe4_expression & !valid_trem2_r47h_variant ~ "Trem2-R47H variant",
      !valid_apoe4_expression & !valid_trem2_r47h_variant ~ "APOE4 expression + Trem2-R47H variant"
    )
  )

plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = is_load1)) +
  geom_point() + ggtitle("Jax.IU.Pitt_APOE4.Trem2.R47H")
print(plt)


## Jax.IU.Pitt_LOAD2.PrimaryScreen ---------------------------------------------

meta_load2 <- subset(metadata_all, study_name == "Jax.IU.Pitt_LOAD2.PrimaryScreen")

# LOAD2 doesn't follow the same nomenclature as LOAD1 genotypes
valid_apoe <- validate_APOE4_KI(meta_load2, counts, symbol_map,
                                genotype_pattern = "LOAD2")
valid_hAbeta <- validate_hAbeta_KI(meta_load2, geno_info,
                                   genotype_pattern = "LOAD2")
valid_trem2 <- validate_Trem2_R47H(meta_load2, geno_info,
                                   genotype_pattern = "LOAD2")

valid_all <- valid_apoe$valid |>
  merge(valid_hAbeta$valid, by = "merged_file_specimenID") |>
  merge(valid_hAbeta$valid, by = "merged_file_specimenID") |>
  mutate(valid = valid.x & valid.y & valid) |>
  select(merged_file_specimenID, valid)

# TODO Il1rap, Il34, Ptprb

# TODO 3 of the 4 WT samples that fail from the LOAD2 study express APOE at ~20
# CPM and should probably pass instead. On the other hand, these 3 look like
# outliers on a PCA

valid_samples <- rbind(valid_samples, valid_all)


## UCI_3xTg-AD -----------------------------------------------------------------

meta_3x <- subset(metadata_all, study_name == "UCI_3xTg-AD")

valid_3x <- validate_3x(meta_3x, geno_info, counts, symbol_map)
valid_samples <- rbind(valid_samples, valid_3x$valid)


## UCI_ABCA7 -------------------------------------------------------------------

meta_abca7 <- subset(metadata_all, study_name == "UCI_ABCA7")

valid_5x <- validate_5x(meta_abca7, geno_info, counts, symbol_map)
valid_abca7 <- validate_Abca7(meta_abca7, geno_info)

valid_all <- merge(valid_5x$valid, valid_abca7$valid, by = "merged_file_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(merged_file_specimenID, valid)

valid_samples <- rbind(valid_samples, valid_all)


## UCI_Clu-h2kbKI --------------------------------------------------------------

meta_clu <- subset(metadata_all, study_name == "UCI_Clu-h2kbKI")

valid_5x <- validate_5x(meta_clu, geno_info, counts, symbol_map)
valid_clu <- validate_CLU_KI(meta_clu, counts, symbol_map)

# There aren't currently any mismatches in this data set
valid_all <- merge(valid_5x$valid, valid_clu$valid, by = "merged_file_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(merged_file_specimenID, valid)

valid_samples <- rbind(valid_samples, valid_all)


## UCI_PrimaryScreen - Abi3-S209F ----------------------------------------------

# Not currently included in analysis

#meta_abi3 <- subset(metadata_all, study_name == "UCI_PrimaryScreen")
#valid_abi3 <- validate_Abi3(meta_abi3, geno_info)

#valid_samples <- valid_samples |>
#  merge(select(valid_abi3, merged_file_specimenID, valid_abi3_variant),
#        all = TRUE)


## UCI_PrimaryScreen - hAbeta-KI -----------------------------------------------

# Not currently included in analysis

#meta_hAbeta <- subset(metadata_all, study_name == "UCI_PrimaryScreen")
#valid_hAbeta <- validate_hAbeta_KI(meta_hAbeta, geno_info)

#valid_samples <- valid_samples |>
#  merge(select(valid_hAbeta, merged_file_specimenID, valid_abeta_ki_variant),
#        all = TRUE)


## UCI_PrimaryScreen - Picalm-H458R --------------------------------------------

# Not currently included in analysis

# Variant for Picalm-H458R on the UCI_PrimaryScreen study seems to be
# unreliable, as the variant is being detected in only one of the Picalm samples
# but deletions at that location are found for multiple samples from unrelated
# studies. Therefore we will not validate based on variant calling. It may be
# difficult to detect that variant because it is 3 consecutive base changes. We
# will re-evaluate if a larger Picalm study gets added.

# Additionally, there is a small difference between genotypes (~2 CPM) but not
# enough to be able to confidently validate by expression, so we currently have
# no way to validate the genotype of these mice.

#valid_samples$valid_Picalm <- TRUE #"Unknown"


## Trem2-KO --------------------------------------------------------------------

# Not currently used in analysis

# Trem2-KO isn't reliably detectable from variant calling due to the length of
# base pair deletion in the mutant mice, so we do not validate based on variant
# calling.

# Trem2-KO mice should not express Trem2. Trem2-R47H samples are removed from
# this comparison so only KO and WT samples are compared.
#studies_trem2ko <- subset(metadata_all, grepl("Trem2-KO", genotype))
#meta_trem2ko <- subset(metadata_all, study_name %in% studies_trem2ko$study_name &
#                         !grepl("R47H", genotype))

#valid_trem2ko <- validate_Trem2_KO(meta_trem2ko, counts, symbol_map)

#valid_samples <- valid_samples |>
#  merge(select(valid_trem2ko, merged_file_specimenID, valid_trem2_ko_expression),
#        all = TRUE)


## UCI_Trem2-R47H_NSS ----------------------------------------------------------

meta_trem2_nss <- subset(metadata_all, study_name == "UCI_Trem2-R47H_NSS")

valid_5x <- validate_5x(meta_trem2_nss, geno_info, counts, symbol_map)
valid_trem2_nss <- validate_Trem2_R47H(meta_trem2_nss, geno_info)

valid_all <- merge(valid_5x$valid, valid_trem2_nss$valid, by = "merged_file_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(merged_file_specimenID, valid)

valid_samples <- rbind(valid_samples, valid_all)


# Combine all validation results -----------------------------------------------

valid_samples <- merge(valid_samples, sex_matches)

stopifnot(nrow(valid_samples) == nrow(metadata_all))
stopifnot(all(valid_samples$merged_file_specimenID %in% metadata_all$merged_file_specimenID))
stopifnot(all(metadata_all$merged_file_specimenID %in% valid_samples$merged_file_specimenID))

valid_samples$validated <- rowSums(valid_samples[, c("valid_sex", "valid")]) == 2

# TODO better output now that I've removed details from the valid_samples df
print(paste(sum(!valid_samples$validated), "samples failed validation:"))
print(subset(valid_samples, !validated))

write.csv(valid_samples, file.path("output", "Model_AD_valid_samples.csv"),
          row.names = FALSE)

# TODO double-check Jax ages match this code
meta_stats <- merge(metadata_all, valid_samples) %>%
  # Fix a few age timepoints for Jax studies
  mutate(
    ageDeath = round(ageDeath),
    ageDeath = case_when(
      study_name == "Jax.IU.Pitt_5XFAD" & ageDeath == 11 ~ 12,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath <= 4 ~ 4,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 4 & ageDeath <= 10 ~ 8,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 10 & ageDeath < 20 ~ 12,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 20 ~ 24,
      .default = ageDeath
    )
  ) %>%
  dplyr::rename(age = ageDeath)

stats <- meta_stats %>%
  group_by(study_name, sex, age, genotype) %>%
  summarize(total_samples = n(),
            valid_samples = sum(validated),
            mismatched_samples = sum(!validated),
            .groups = "drop")

write.csv(stats, "study_validated_samples_stats.csv", row.names = FALSE)

# TODO put on Synapse
