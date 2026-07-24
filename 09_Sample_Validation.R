library(synapser)
library(synapserutils)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

source("validation_functions.R")
source("util_functions.R")

# Set up -----------------------------------------------------------------------

file_syn_ids <- config::get("file_syn_ids")
studies <- config::get("studies")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url"),
                 "/blob/main/09_Sample_Validation.R")
tmp_dir <- file.path("output", "tmp")


# Utility functions ------------------------------------------------------------

# Get a data frame with specific genes, merged with metadata
make_counts_df <- function(metadata, counts, symbol_map, genes) {
  genes_sub <- subset(symbol_map, gene_symbol %in% genes)
  rownames(genes_sub) <- genes_sub$ensembl_gene_id

  counts_sub <- counts[genes_sub$ensembl_gene_id, ]
  rownames(counts_sub) <- genes_sub[rownames(counts_sub), "gene_symbol"]

  counts_df <- as.data.frame(t(counts_sub)) |>
    merge(metadata, by.x = "row.names", by.y = "unique_specimenID") |>
    dplyr::rename(unique_specimenID = Row.names)
}


score_quality <- function(quality) {
  case_when(
    is.na(quality) ~ "deletion",
    quality < 20 ~ "low",
    quality >= 20 & quality < 100 ~ "moderate",
    quality >= 100 ~ "high"
  )
}


get_variant_mismatches <- function(metadata, geno_info, genotype_pattern,
                                   mutation_match = NULL,
                                   gene_symbol_match = NULL,
                                   total_positions = 0,
                                   rm_na = TRUE) {
  carrier_samples <- subset(metadata, grepl(genotype_pattern, genotype))

  # Get every sample in the same studies as the carrier samples
  study_samples <- subset(metadata, study %in% carrier_samples$study)

  geno_sub <- subset(geno_info,
                     unique_specimenID %in% study_samples$unique_specimenID &
                       (gene_symbol %in% gene_symbol_match |
                          mutation %in% mutation_match)) |>
    group_by(study, unique_specimenID, specimenID, genotype, gene_symbol) |>
    summarize(count = sum(evidence != "low"),
              avg_quality = mean(quality, na.rm = TRUE),
              .groups = "drop") |>
    mutate(evidence = score_quality(avg_quality),
           is_carrier = grepl(genotype_pattern, genotype))

  # Add in samples that don't show up in geno_sub
  missing_samples <- study_samples |>
    subset(!(unique_specimenID %in% geno_sub$unique_specimenID)) |>
    group_by(study, unique_specimenID, specimenID, genotype) |>
    reframe(gene_symbol = unique(geno_sub$gene_symbol),
            count = 0, avg_quality = 0, evidence = "none",
            is_carrier = grepl(genotype_pattern, genotype))

  geno_sub <- rbind(geno_sub, missing_samples)

  geno_final <- geno_sub |>
    group_by(study, unique_specimenID, specimenID, genotype, is_carrier) |>
    summarize(total_found = sum(count),
              evidence = paste(unique(evidence), collapse = ","),
              .groups = "drop")

  if (rm_na) {
    geno_final$total_found[geno_final$evidence == "deletion"] <- 0
  }

  geno_final <- geno_final |>
    mutate(est_genotype =
             ifelse(total_found >= total_positions, "carrier",
                    ifelse(total_found == 0 & evidence == "none", "noncarrier",
                           "ambiguous")))

  return(geno_final)
}


reformat_details <- function(details_df, symbol_map) {
  details_df <- details_df |>
    select(study, unique_specimenID, genotype,
           any_of(c("total_found", symbol_map$gene_symbol)),
           matches("valid_.*_variant"), matches("valid_.*_expression"))

  if ("total_found" %in% colnames(details_df)) {
    details_df <- dplyr::rename(details_df, total_variants_found = total_found)
  }

  return(details_df)
}


# Load counts and metadata -----------------------------------------------------

meta_list <- get_all_metadata()
meta_list <- meta_list[studies] |>
  lapply("[[", "data")

# Combine into one data frame and bin the continuous ageDeath values in Jax studies
metadata_all <- do.call(rbind, meta_list) |>
  bin_jax_ages()

symbol_map_file <- synGet(file_syn_ids$symbol_map,
                          downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

symbol_map <- read.csv(symbol_map_file$path) |>
  arrange(ensembl_gene_id)

# Read in all counts files
cat("Reading counts files...\n")
counts_list <- get_all_counts_files(studies,
                                    meta_list,
                                    symbol_map,
                                    count_type = "gene_counts") |>
  lapply("[[", "counts")

counts <- do.call(cbind, counts_list)

# Convert counts to CPM so samples are comparable to each other
lib_sizes <- colSums(counts)
counts <- sweep(counts, 2, lib_sizes, "/") * 1e6

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, unique_specimenID %in% colnames(counts))
counts <- counts[, metadata_all$unique_specimenID]


# Load intervals file ----------------------------------------------------------

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

intervals <- read.delim(intervals_file$path, header = FALSE) |>
  dplyr::rename(chromosome = V1, start = V2, end = V3, mutation = V4) |>
  rowwise() |>
  # adds as many new rows as necessary to cover the full range of positions
  reframe(chromosome = chromosome,
          mutation = mutation,
          position = (start+1):(end-1))

# Load VCF files ---------------------------------------------------------------

# Find the combined VCF files on Synapse
vcf_files <- syn_get_unique_children("genotype_validation")
vcf_files <- vcf_files[grepl("combined_vcf", names(vcf_files))]

# Merge with intervals data to get which gene belongs to each position
geno_info <- lapply(vcf_files, function(study_file) {
  syn_file <- synGet(study_file$id,
                     downloadLocation = tmp_dir,
                     ifcollision = "overwrite.local")

  data <- read.csv(syn_file$path) |>
    merge(intervals) |>
    mutate(gene_symbol = str_replace(mutation, "[-\\*].*", ""))

  return(data)
})

geno_info <- do.call(rbind, geno_info)
geno_info$quality[geno_info$quality == "."] <- NA
geno_info$quality <- as.numeric(geno_info$quality)
geno_info <- merge(geno_info, metadata_all) |>
  mutate(evidence = score_quality(quality))


# Set up the valid samples list ------------------------------------------------

valid_samples_list <- vector("list", length(unique(metadata_all$study)))
names(valid_samples_list) <- unique(metadata_all$study)


# Validation of mouse sex ------------------------------------------------------

xy_df <- make_counts_df(metadata_all, counts, symbol_map,
                        c("Xist", "Eif2s3y", "Ddx3y", "Kdm5d")) |>
  # Use a small threshold for Xist for females -- there is one sample that
  # expresses near-zero counts of both Y-genes and Xist.
  mutate(
    # Take the geometric mean (log-scale), including a pseudocount
    mean_y = rowMeans(log(cbind(Eif2s3y, Ddx3y, Kdm5d) + 0.5)),
    mean_y = exp(mean_y) - 0.5,
    # Expression of Y-related genes should be zero for females, so we use
    # anything > 1 CPM for males and assume < 1 CPM might be noise.
    est_male = mean_y >= 1,
    est_female = Xist > 10 & mean_y < 1
  )

# Mark valid/invalid conditions
xy_df$valid_sex <- FALSE
xy_df$valid_sex[xy_df$sex == "female" &
                  xy_df$est_female == TRUE &
                  xy_df$est_male == FALSE] <- TRUE
xy_df$valid_sex[xy_df$sex == "male" &
                  xy_df$est_female == FALSE &
                  xy_df$est_male == TRUE] <- TRUE

# Linear scale
ggplot(xy_df, aes(x = Xist, y = mean_y, color = sex)) +
  geom_point() +
  geom_hline(yintercept = 1, linewidth = 0.5, color = "orange") +
  theme_bw() +
  facet_wrap(~study)

# Log scale - Adding pseudocount to avoid inf values
ggplot(xy_df, aes(x = log2(Xist+0.5), y = log2(mean_y+0.5), color = sex)) +
  geom_point() +
  geom_hline(yintercept = log2(1.5), linewidth = 0.5, color = "orange") +
  theme_bw() +
  facet_wrap(~study)

# Mark valid samples in the data frame
sex_matches <- select(xy_df, unique_specimenID, valid_sex)

# For printing
mismatches <- subset(xy_df, !valid_sex) |>
  select(study, individualID, specimenID, sex, est_female, est_male,
         Xist, Eif2s3y, Ddx3y, Kdm5d, mean_y) |>
  dplyr::rename(reported_sex = sex) |>
  mutate(across(c(Xist, Ddx3y, Eif2s3y, Kdm5d, mean_y), ~ round(.x, 2)))

mismatches$est_sex <- "unknown"
mismatches$est_sex[mismatches$est_female & !mismatches$est_male] <- "female"
mismatches$est_sex[mismatches$est_male & !mismatches$est_female] <- "male"

print(paste(nrow(mismatches), "samples have mismatched sex:"))
print(select(mismatches, study, specimenID, reported_sex, est_sex,
             Xist, Eif2s3y, Ddx3y, Kdm5d, mean_y))

# The two Jax.IU.Pitt_LOAD2.PrimaryScreen samples that fail (both female) have a
# mean_y of 3.64 and 11.06, which is much lower than male samples but higher
# than all other female samples. These samples also express a small amount of
# APOE (~20 CPM) but are supposed to be APOE_WT. Combined, this suggests
# contamination of these samples.
#
# The four Jax.IU.Pitt_LOAD2 samples that fail (all female) have mean_y ranging
# from 1.42 to 8.36 but Xist > 1200. Possibly two of these could pass with < 3
# CPM mean_y but it's unclear. These samples don't fail expression or variant
# checks.

# Validation of mouse genotype -------------------------------------------------

## Jax.IU.Pitt_5XFAD and UCI_5XFAD ---------------------------------------------

meta_5x <- subset(metadata_all,
                  study %in% c("Jax.IU.Pitt_5XFAD", "UCI_5XFAD"))

valid_5x <- validate_5x(meta_5x, geno_info, counts, symbol_map)

# TODO temporarily duplicated data
valid_samples_list[["Jax.IU.Pitt_5XFAD"]] <- valid_5x$valid
valid_samples_list[["UCI_5XFAD"]] <- valid_5x$valid

# Tmp PCA

counts_5x <- log2(counts[, meta_5x$unique_specimenID] + 0.5)

for (study_name in unique(meta_5x$study)) {
  meta_study <- subset(meta_5x, study == study_name) |>
    merge(valid_5x$detail)
  vars <- matrixStats::rowVars(as.matrix(counts_5x[, meta_study$unique_specimenID]))
  hv <- names(sort(vars, decreasing = TRUE))[1:2000]
  pc <- prcomp(t(counts_5x[hv, meta_study$unique_specimenID]))
  pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                   by.y = "unique_specimenID") |>
    dplyr::rename(unique_specimenID = Row.names) |>
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
    geom_point() + ggtitle(study_name)
  print(plt)
}

# NOTE: GT19_12887 (non-carrier) from Jax 5X study has all 6 variants detected
# but expresses only a small amount of APP (11.28 CPM) and PSEN1 (0.85 CPM),
# and clusters with other non-carriers in a PCA


## Jax.IU.Pitt_APOE4.Trem2.R47H ------------------------------------------------

# APOE4 KI can be validated with gene expression only, and Trem2-R47H can be
# validated with variant calling.
meta_load1 <- subset(metadata_all, study == "Jax.IU.Pitt_APOE4.Trem2.R47H")

valid_apoe <- validate_APOE4_KI(meta_load1, counts, symbol_map)
valid_trem2 <- validate_Trem2_R47H(meta_load1, geno_info)

valid_all <- merge(valid_apoe$valid, valid_trem2$valid,
                   by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(unique_specimenID, valid)

valid_samples_list[["Jax.IU.Pitt_APOE4.Trem2.R47H"]] <- valid_all

# Tmp PCA

meta_study <- subset(meta_load1, !is.na(genotype)) |>
  merge(valid_apoe$detail) |>
  merge(valid_trem2$detail)

counts_apoe <- log2(counts[, meta_study$unique_specimenID] + 0.5)

vars <- matrixStats::rowVars(as.matrix(counts_apoe[, meta_study$unique_specimenID]))
hv <- names(sort(vars, decreasing = TRUE))[1:2000]

pc <- prcomp(t(counts_apoe[hv, meta_study$unique_specimenID]))

pc_plot <- merge(pc$x[, paste0("PC", 1:10)], meta_study, by.x = "row.names",
                 by.y = "unique_specimenID") |>
  dplyr::rename(unique_specimenID = Row.names) |>
  mutate(
    is_load1 = grepl("APOE4-KI_(homo|hetero)", genotype),
    mismatch_type = case_when(
      valid_apoe4_expression & valid_trem2_r47h_variant ~ "none",
      !valid_apoe4_expression & valid_trem2_r47h_variant ~ "APOE4 expression",
      valid_apoe4_expression & !valid_trem2_r47h_variant ~ "Trem2-R47H variant",
      !valid_apoe4_expression & !valid_trem2_r47h_variant ~ "APOE4 expression + Trem2-R47H variant"
    )
  )

# PC1 separates by sequencingBatch, PC2 separates by sex
plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = is_load1)) +
  geom_point() + ggtitle("Jax.IU.Pitt_APOE4.Trem2.R47H")
print(plt)

plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = specimenID %in% mismatches$specimenID, shape = sex)) +
  geom_point() + ggtitle("Jax.IU.Pitt_APOE4.Trem2.R47H")
print(plt)

# PC3 separates by genotype
plt <- ggplot(pc_plot, aes(x = PC2, y = PC3, color = mismatch_type, shape = is_load1)) +
  geom_point() + ggtitle("Jax.IU.Pitt_APOE4.Trem2.R47H")
print(plt)

# PC4 somewhat separates by ageDeath
plt <- ggplot(pc_plot, aes(x = PC2, y = PC3, color = factor(ageDeath), shape = is_load1)) +
  geom_point() + ggtitle("Jax.IU.Pitt_APOE4.Trem2.R47H")
print(plt)

# Expression plots

counts_load1 <- make_counts_df(meta_load1, counts, symbol_map,
                               c("APOE", "Apoe", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load1, aes(x = genotype, y = expr, color = valid, shape = genotype)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")


## Jax.IU.Pitt_LOAD2 -----------------------------------------------------------

# APOE4 KI can be validated with gene expression only, and Trem2-R47H and APP-KI
# can be validated with variant calling.
meta_load2 <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2")

# TODO there is some ambiguous WT APOE expression. All 3 ambiguous samples
# express Apoe and Trem2 at a rate comparable to other APOE non-carriers and should
# maybe pass instead of failing?
valid_apoe <- validate_APOE4_KI(meta_load2, counts, symbol_map)
valid_trem2 <- validate_Trem2_R47H(meta_load2, geno_info)
valid_hAbeta <- validate_hAbeta_KI(meta_load2, geno_info,
                                   genotype_pattern = "hAPP-3/3_homo|hetero")

valid_all <- merge(valid_apoe$valid, valid_trem2$valid,
                   by = "unique_specimenID") |>
  merge(valid_hAbeta$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y & valid) |>
  select(unique_specimenID, valid)

valid_samples_list[["Jax.IU.Pitt_LOAD2"]] <- valid_all

# Tmp PCA

meta_study <- subset(meta_load2, !is.na(genotype)) |>
  merge(select(valid_apoe$detail, unique_specimenID, valid_apoe4_expression)) |>
  merge(select(valid_trem2$detail, unique_specimenID, valid_trem2_r47h_variant)) |>
  merge(select(valid_hAbeta$detail, unique_specimenID, valid_hAbeta_variant))

counts_load2 <- log2(counts[, meta_study$unique_specimenID] + 0.5)

vars <- matrixStats::rowVars(as.matrix(counts_load2[, meta_study$unique_specimenID]))
hv <- names(sort(vars, decreasing = TRUE))[1:2000]

pc <- prcomp(t(counts_load2[hv, meta_study$unique_specimenID]))

pc_plot <- merge(pc$x[, paste0("PC", 1:10)], meta_study, by.x = "row.names",
                 by.y = "unique_specimenID") |>
  dplyr::rename(unique_specimenID = Row.names) |>
  mutate(
    has_apoe4 = grepl("APOE4-KI_(homo|hetero)", genotype),
    mismatch_type = case_when(
      valid_apoe4_expression & valid_trem2_r47h_variant & valid_hAbeta_variant ~ "none",
      !valid_apoe4_expression & valid_trem2_r47h_variant & valid_hAbeta_variant ~ "APOE4 expression",
      valid_apoe4_expression & !valid_trem2_r47h_variant & valid_hAbeta_variant ~ "Trem2-R47H variant",
      valid_apoe4_expression & valid_trem2_r47h_variant & !valid_hAbeta_variant ~ "hAPP variant",
      .default = "Multiple mismatches"
    )
  )

# PC1 separates by sequencingBatch and somewhat by ageDeath
plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = has_apoe4)) +
  geom_point() + ggtitle("Jax.IU.Pitt_LOAD2")
print(plt)

# PC3 separates by sex -- highlight the sex mismatches
plt <- ggplot(pc_plot, aes(x = PC2, y = PC3, color = specimenID %in% mismatches$specimenID, shape = sex)) +
  geom_point() + ggtitle("Jax.IU.Pitt_LOAD2")
print(plt)

# PC7 separates by genotype
plt <- ggplot(pc_plot, aes(x = PC1, y = PC7, color = mismatch_type, shape = has_apoe4)) +
  geom_point() + ggtitle("Jax.IU.Pitt_LOAD2")
print(plt)

# Expression plots

counts_load2 <- make_counts_df(meta_load2, counts, symbol_map,
                               c("APOE", "Apoe", "APP", "App", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, APP, App, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load2, aes(x = genotype, y = expr, color = valid, shape = genotype)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")


## Jax.IU.Pitt_LOAD2.PrimaryScreen ---------------------------------------------

meta_load2pri <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2.PrimaryScreen")

# LOAD2 doesn't follow the same nomenclature as LOAD1 genotypes
valid_apoe <- validate_APOE4_KI(meta_load2pri, counts, symbol_map,
                                genotype_pattern = "LOAD2")
valid_hAbeta <- validate_hAbeta_KI(meta_load2pri, geno_info,
                                   genotype_pattern = "LOAD2")
valid_trem2 <- validate_Trem2_R47H(meta_load2pri, geno_info,
                                   genotype_pattern = "LOAD2")

valid_all <- valid_apoe$valid |>
  merge(valid_hAbeta$valid, by = "unique_specimenID") |>
  merge(valid_trem2$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y & valid) |>
  select(unique_specimenID, valid)

# TODO Il1rap, Il34, Ptprb

valid_samples_list[["Jax.IU.Pitt_LOAD2.PrimaryScreen"]] <- valid_all

# Tmp PCA

meta_study <- subset(meta_load2pri, !is.na(genotype)) |>
  merge(valid_apoe$detail) |>
  merge(valid_hAbeta$detail) |>
  merge(valid_trem2$detail,
        by = c("study", "specimenID", "unique_specimenID", "genotype"))

counts_load2 <- log2(counts[, meta_study$unique_specimenID] + 0.5)

vars <- matrixStats::rowVars(as.matrix(counts_load2))
hv <- names(sort(vars, decreasing = TRUE))[1:2000]

pc <- prcomp(t(counts_load2[hv, ]))

pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                 by.y = "unique_specimenID") |>
  dplyr::rename(unique_specimenID = Row.names) |>
  mutate(
    is_load2 = grepl("LOAD2", genotype),
    mismatch_type = case_when(
      # There aren't any hAbeta mismatches so it isn't accounted for here
      valid_apoe4_expression & valid_trem2_r47h_variant ~ "none",
      !valid_apoe4_expression & valid_trem2_r47h_variant ~ "APOE4 expression",
      valid_apoe4_expression & !valid_trem2_r47h_variant ~ "Trem2-R47H variant",
      !valid_apoe4_expression & !valid_trem2_r47h_variant ~ "APOE4 expression + Trem2-R47H variant"
    )
  )

plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = is_load2)) +
  geom_point() + ggtitle("Jax.IU.Pitt_LOAD2.PrimaryScreen")
print(plt)

counts_load2 <- make_counts_df(meta_load2pri, counts, symbol_map,
                               c("APOE", "Apoe", "APP", "App", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, APP, App, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load2, aes(x = genotype, y = expr, color = valid)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")

# TODO The three non-carrier samples that express a small amount of APOE express
# Apoe and Trem2 at levels comparable to other non-carriers. Should they pass
# instead?


## UCI_3xTg-AD -----------------------------------------------------------------

meta_3x <- subset(metadata_all, study == "UCI_3xTg-AD")

valid_3x <- validate_3x(meta_3x, geno_info, counts, symbol_map)
valid_samples_list[["UCI_3xTg-AD"]] <- valid_3x$valid


## UCI_ABCA7 -------------------------------------------------------------------

meta_abca7 <- subset(metadata_all, study == "UCI_ABCA7")

valid_5x <- validate_5x(meta_abca7, geno_info, counts, symbol_map)
valid_abca7 <- validate_Abca7(meta_abca7, geno_info)

valid_all <- merge(valid_5x$valid, valid_abca7$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(unique_specimenID, valid)

valid_samples_list[["UCI_ABCA7"]] <- valid_all

# Sample 11451lh has all 6 detected variants but only expresses 11.2 CPM APP and
# 0.5 CPM PSEN1. Its APP expression is an outlier compared to other
# 5XFAD_noncarriers but not nearly high enough to suggest the actual 5XFAD
# genotype. To be cautious, this sample is dropped.

counts_abca7 <- make_counts_df(meta_abca7, counts, symbol_map,
                               c("APP", "App", "PSEN1", "Psen1", "Abca7")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, Abca7), names_to = "gene", values_to = "expr")

ggplot(counts_abca7, aes(x = genotype, y = expr, color = valid, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


## UCI_Clu-h2kbKI --------------------------------------------------------------

meta_clu <- subset(metadata_all, study == "UCI_Clu-h2kbKI")

valid_5x <- validate_5x(meta_clu, geno_info, counts, symbol_map)
valid_clu <- validate_CLU_KI(meta_clu, counts, symbol_map)

# No mismatches in this data set
valid_all <- merge(valid_5x$valid, valid_clu$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(unique_specimenID, valid)

valid_samples_list[["UCI_Clu-h2kbKI"]] <- valid_all

counts_clu <- make_counts_df(meta_clu, counts, symbol_map,
                             c("APP", "App", "PSEN1", "Psen1", "CLU", "Clu")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, CLU, Clu), names_to = "gene", values_to = "expr")

ggplot(counts_clu, aes(x = genotype, y = expr, color = genotype, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


## UCI_Bin1K358R ---------------------------------------------------------------

meta_bin1 <- subset(metadata_all, study == "UCI_Bin1K358R")

valid_5x <- validate_5x(meta_bin1, geno_info, counts, symbol_map)
valid_bin1 <- validate_Bin1(meta_bin1, geno_info, counts, symbol_map)

# No mismatches in this data set
valid_all <- merge(valid_5x$valid, valid_bin1$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(unique_specimenID, valid)

valid_samples_list[["UCI_Bin1K358R"]] <- valid_all

counts_bin1 <- make_counts_df(meta_bin1, counts, symbol_map,
                              c("APP", "App", "PSEN1", "Psen1", "Bin1")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, Bin1), names_to = "gene", values_to = "expr")

ggplot(counts_bin1, aes(x = genotype, y = expr, color = genotype, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


## UCI_PrimaryScreen - Abi3-S209F ----------------------------------------------

# Not currently included in analysis

#meta_abi3 <- subset(metadata_all, study == "UCI_PrimaryScreen")
#valid_abi3 <- validate_Abi3(meta_abi3, geno_info)


## UCI_PrimaryScreen - hAbeta-KI -----------------------------------------------

# Not currently included in analysis

#meta_hAbeta <- subset(metadata_all, study == "UCI_PrimaryScreen")
#valid_hAbeta <- validate_hAbeta_KI(meta_hAbeta, geno_info)


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
#meta_trem2ko <- subset(metadata_all, study %in% studies_trem2ko$study &
#                         !grepl("R47H", genotype))

#valid_trem2ko <- validate_Trem2_KO(meta_trem2ko, counts, symbol_map)


## UCI_Trem2-R47H_NSS ----------------------------------------------------------

meta_trem2_nss <- subset(metadata_all, study == "UCI_Trem2-R47H_NSS")

valid_5x <- validate_5x(meta_trem2_nss, geno_info, counts, symbol_map)
valid_trem2_nss <- validate_Trem2_R47H(meta_trem2_nss, geno_info)

valid_all <- merge(valid_5x$valid, valid_trem2_nss$valid, by = "unique_specimenID") |>
  mutate(valid = valid.x & valid.y) |>
  select(unique_specimenID, valid)

valid_samples_list[["UCI_Trem2-R47H_NSS"]] <- valid_all

counts_trem2 <- make_counts_df(meta_trem2_nss, counts, symbol_map,
                               c("APP", "App", "PSEN1", "Psen1", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_trem2, aes(x = genotype, y = expr, color = valid, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


# Combine all validation results -----------------------------------------------

valid_samples <- do.call(rbind, valid_samples_list) |>
  merge(sex_matches) |>
  distinct()

stopifnot(nrow(valid_samples) == nrow(metadata_all))
stopifnot(all(valid_samples$unique_specimenID %in% metadata_all$unique_specimenID))
stopifnot(all(metadata_all$unique_specimenID %in% valid_samples$unique_specimenID))

valid_samples$validated <- rowSums(valid_samples[, c("valid_sex", "valid")]) == 2

# TODO better output now that I've removed details from the valid_samples df
print(paste(sum(!valid_samples$validated), "samples failed validation:"))
print(subset(valid_samples, !validated))

# TODO maybe split "valid" into "valid_genotype" and "valid_expression" for
# readability
# TODO rename to "genotype_validated_samples.csv" and put in Metadata folder?
write.csv(valid_samples, file.path("output", "Model_AD_valid_samples.csv"),
          row.names = FALSE, quote = FALSE)

meta_stats <- merge(metadata_all, valid_samples) |>
  dplyr::rename(age = ageDeath)

stats <- meta_stats |>
  group_by(study, sex, age, genotype) |>
  summarize(total_samples = n(),
            valid_samples = sum(validated),
            mismatched_samples = sum(!validated),
            .groups = "drop")

write.csv(stats, "study_validated_samples_stats.csv",
          row.names = FALSE, quote = FALSE)

# TODO put on Synapse
