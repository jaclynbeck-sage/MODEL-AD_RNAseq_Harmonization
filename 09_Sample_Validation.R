library(synapser)
library(synapserutils)
library(dplyr)
library(stringr)
library(purrr)

source(file.path("functions", "validation_functions.R"))
source(file.path("functions", "plotting_functions.R"))
source(file.path("functions", "util_functions.R"))

# Set up -----------------------------------------------------------------------

file_syn_ids <- config::get("file_syn_ids")
studies <- config::get("studies")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url"),
                 "/blob/main/09_Sample_Validation.R")

tmp_dir <- file.path("output", "tmp")
dir.create(tmp_dir, showWarnings = FALSE)


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
geno_list <- lapply(vcf_files, function(study_file) {
  syn_file <- synGet(study_file$id,
                     downloadLocation = tmp_dir,
                     ifcollision = "overwrite.local")

  data <- read.csv(syn_file$path) |>
    merge(intervals) |>
    mutate(gene_symbol = str_replace(mutation, "[-\\*].*", ""),
           # Quality for non-deletions is numeric data, while quality for
           # deletions is ".". This turns "." into NA for easier manipulation
           # later.
           quality = suppressWarnings(as.numeric(quality)),
           evidence = score_quality(quality))

  return(data)
})

# Merge in unique_specimenIDs from metadata
geno_info <- do.call(rbind, geno_list) |>
  merge(select(metadata_all, study, unique_specimenID, specimenID))


# Set up the valid samples list ------------------------------------------------

valid_samples_list <- vector("list", length(unique(metadata_all$study)))
names(valid_samples_list) <- unique(metadata_all$study)


# Validation of mouse sex ------------------------------------------------------

sex_matches <- validate_sex(metadata_all, counts, symbol_map)

# Linear scale
ggplot(sex_matches, aes(x = Xist, y = mean_y, color = sex)) +
  geom_point() +
  geom_hline(yintercept = 1, linewidth = 0.5, color = "orange") +
  theme_bw() +
  facet_wrap(~study)

# Log scale - Adding pseudocount to avoid inf values
ggplot(sex_matches, aes(x = log2(Xist+0.5), y = log2(mean_y+0.5), color = sex)) +
  geom_point() +
  geom_hline(yintercept = log2(1.5), linewidth = 0.5, color = "orange") +
  theme_bw() +
  facet_wrap(~study)

# NOTE:
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

cat("*** Jax.IU.Pitt_5XFAD, UCI_5XFAD ***", "\n", "\n")
meta_5x <- subset(metadata_all,
                  study %in% c("Jax.IU.Pitt_5XFAD", "UCI_5XFAD"))

valid_5x <- validate_5x(meta_5x, geno_info, counts, symbol_map)

valid_samples_list[["Jax.IU.Pitt_5XFAD"]] <- subset(valid_5x, study == "Jax.IU.Pitt_5XFAD")
valid_samples_list[["UCI_5XFAD"]] <- subset(valid_5x, study == "UCI_5XFAD")

# PCA
for (study_name in unique(meta_5x$study)) {
  pc_plot <- calculate_pca(subset(valid_5x, study == study_name), counts,
                           "5XFAD_carrier") |>
    dplyr::rename(has_5x = has_mutation)

  plot_pc_grid(pc_plot, study_name, "has_5x")
}

# NOTE: GT19_12887 (non-carrier) from Jax 5X study has all 6 variants detected
# but expresses only a small amount of APP (11.28 CPM) and PSEN1 (0.85 CPM),
# and clusters with other non-carriers in a PCA


## Jax.IU.Pitt_APOE4.Trem2.R47H ------------------------------------------------

cat("*** Jax.IU.Pitt_APOE4.Trem2.R47H ***", "\n", "\n")

# APOE4 KI can be validated with gene expression only, and Trem2-R47H can be
# validated with variant calling.
meta_load1 <- subset(metadata_all, study == "Jax.IU.Pitt_APOE4.Trem2.R47H")

valid_apoe <- validate_APOE4_KI(meta_load1, counts, symbol_map)
valid_trem2 <- validate_Trem2_R47H(meta_load1, geno_info)

valid_all <- merge_validation_dfs(valid_apoe, valid_trem2)

valid_samples_list[["Jax.IU.Pitt_APOE4.Trem2.R47H"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "APOE4-KI_(homo|hetero)") |>
  rename(has_apoe = has_mutation)

plot_pc_grid(pc_plot, "Jax.IU.Pitt_APOE4.Trem2.R47H", "has_apoe")

# Expression plots

counts_load1 <- make_counts_df(meta_load1, counts, symbol_map,
                               c("APOE", "Apoe", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load1, aes(x = genotype, y = expr, color = valid_expression, shape = genotype)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")


## Jax.IU.Pitt_LOAD2 -----------------------------------------------------------

cat("*** Jax.IU.Pitt_LOAD2 ***", "\n", "\n")

# APOE4 KI can be validated with gene expression only, and Trem2-R47H and APP-KI
# can be validated with variant calling.
meta_load2 <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2")

# TODO there is some ambiguous WT APOE expression. All 3 ambiguous samples
# express Apoe and Trem2 at a rate comparable to other APOE non-carriers and should
# maybe pass instead of failing?
valid_apoe <- validate_APOE4_KI(meta_load2, counts, symbol_map)
valid_trem2 <- validate_Trem2_R47H(meta_load2, geno_info)
valid_hAPP <- validate_hAPP_KI(meta_load2, geno_info)

valid_all <- merge_validation_dfs(valid_apoe, valid_trem2) |>
  merge_validation_dfs(valid_hAPP)

valid_samples_list[["Jax.IU.Pitt_LOAD2"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "APOE4-KI_(homo|hetero)") |>
  dplyr::rename(has_apoe4 = has_mutation)

plot_pc_grid(pc_plot, "Jax.IU.Pitt_LOAD2", "has_apoe4")

# Expression plots

counts_load2 <- make_counts_df(meta_load2, counts, symbol_map,
                               c("APOE", "Apoe", "APP", "App", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, APP, App, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load2, aes(x = genotype, y = expr, color = valid_expression, shape = genotype)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")


## Jax.IU.Pitt_LOAD2.PrimaryScreen ---------------------------------------------

cat("*** Jax.IU.Pitt_LOAD2.PrimaryScreen ***", "\n", "\n")
meta_load2pri <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2.PrimaryScreen")

# LOAD2 doesn't follow the same nomenclature as LOAD1 genotypes
valid_apoe <- validate_APOE4_KI(meta_load2pri, counts, symbol_map,
                                genotype_pattern = "LOAD2")
valid_hAPP <- validate_hAPP_KI(meta_load2pri, geno_info,
                               genotype_pattern = "LOAD2")
valid_trem2 <- validate_Trem2_R47H(meta_load2pri, geno_info,
                                   genotype_pattern = "LOAD2")

valid_all <- merge_validation_dfs(valid_apoe, valid_trem2) |>
  merge_validation_dfs(valid_hAPP)

# TODO Il1rap, Il34, Ptprb

valid_samples_list[["Jax.IU.Pitt_LOAD2.PrimaryScreen"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "LOAD2") |>
  dplyr::rename(is_load2 = has_mutation)

plot_pc_grid(pc_plot, "Jax.IU.Pitt_LOAD2.PrimaryScreen", "is_load2")

# Counts plots

counts_load2 <- make_counts_df(meta_load2pri, counts, symbol_map,
                               c("APOE", "Apoe", "APP", "App", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APOE, Apoe, APP, App, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_load2, aes(x = genotype, y = expr, color = valid_expression)) +
  geom_jitter() +
  facet_grid(rows = vars(gene), cols = vars(ageDeath), scales = "free")

# TODO The three non-carrier samples that express a small amount of APOE express
# Apoe and Trem2 at levels comparable to other non-carriers. Should they pass
# instead?


## UCI_3xTg-AD -----------------------------------------------------------------

cat("*** UCI_3xTg-AD ***", "\n", "\n")
meta_3x <- subset(metadata_all, study == "UCI_3xTg-AD")

valid_3x <- validate_3x(meta_3x, geno_info, counts, symbol_map)
valid_samples_list[["UCI_3xTg-AD"]] <- valid_3x

pc_plot <- calculate_pca(valid_3x, counts, "3xTg-AD_carrier")
plot_pc_grid(pc_plot, "UCI_3xTg-AD")


## UCI_ABCA7 -------------------------------------------------------------------

cat("*** UCI_ABCA7 ***", "\n", "\n")
meta_abca7 <- subset(metadata_all, study == "UCI_ABCA7")

valid_5x <- validate_5x(meta_abca7, geno_info, counts, symbol_map)
valid_abca7 <- validate_Abca7(meta_abca7, geno_info)

valid_all <- merge_validation_dfs(valid_5x, valid_abca7)
valid_samples_list[["UCI_ABCA7"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "Abca7.*_homozygous") |>
  mutate(has_5x = grepl("5XFAD_carrier", genotype))
plot_pc_grid(pc_plot, "UCI_ABCA7", "has_5x")

# Sample 11451lh has all 6 detected variants but only expresses 11.2 CPM APP and
# 0.5 CPM PSEN1. Its APP expression is an outlier compared to other
# 5XFAD_noncarriers but not nearly high enough to suggest the actual 5XFAD
# genotype. To be cautious, this sample is dropped.

counts_abca7 <- make_counts_df(meta_abca7, counts, symbol_map,
                               c("APP", "App", "PSEN1", "Psen1", "Abca7")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, Abca7), names_to = "gene", values_to = "expr")

ggplot(counts_abca7, aes(x = genotype, y = expr, color = valid_expression, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


## UCI_Clu-h2kbKI --------------------------------------------------------------

cat("*** UCI_Clu-h2kbKI ***", "\n", "\n")
meta_clu <- subset(metadata_all, study == "UCI_Clu-h2kbKI")

valid_5x <- validate_5x(meta_clu, geno_info, counts, symbol_map)
valid_clu <- validate_CLU_KI(meta_clu, counts, symbol_map)

# No mismatches in this data set
valid_all <- merge_validation_dfs(valid_5x, valid_clu)
valid_samples_list[["UCI_Clu-h2kbKI"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "Clu.*_homozygous") |>
  mutate(has_5x = grepl("5XFAD_carrier", genotype))
plot_pc_grid(pc_plot, "UCI_Clu-h2kbKI", "has_5x")

# Counts

counts_clu <- make_counts_df(meta_clu, counts, symbol_map,
                             c("APP", "App", "PSEN1", "Psen1", "CLU", "Clu")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, CLU, Clu), names_to = "gene", values_to = "expr")

ggplot(counts_clu, aes(x = genotype, y = expr, color = genotype, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


## UCI_Bin1K358R ---------------------------------------------------------------

cat("*** UCI_Bin1K358R ***", "\n", "\n")
meta_bin1 <- subset(metadata_all, study == "UCI_Bin1K358R")

valid_5x <- validate_5x(meta_bin1, geno_info, counts, symbol_map)
valid_bin1 <- validate_Bin1(meta_bin1, geno_info, counts, symbol_map)

# No mismatches in this data set
valid_all <- merge_validation_dfs(valid_5x, valid_bin1)
valid_samples_list[["UCI_Bin1K358R"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "Bin1.*_homozygous") |>
  mutate(has_5x = grepl("5XFAD_carrier", genotype))
plot_pc_grid(pc_plot, "UCI_Bin1K358R", "has_5x")

# Counts

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
#valid_hAPP <- validate_hAPP(meta_hAbeta, geno_info,
#                            genotype_pattern = "hAbeta-KI_LoxP_homozygous")


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

cat("*** UCI_Trem2-R47H_NSS ***", "\n", "\n")
meta_trem2_nss <- subset(metadata_all, study == "UCI_Trem2-R47H_NSS")

valid_5x <- validate_5x(meta_trem2_nss, geno_info, counts, symbol_map)
valid_trem2_nss <- validate_Trem2_R47H(meta_trem2_nss, geno_info)

valid_all <- merge_validation_dfs(valid_5x, valid_trem2_nss)
valid_samples_list[["UCI_Trem2-R47H_NSS"]] <- valid_all

# PCA

pc_plot <- calculate_pca(valid_all, counts, "Trem2.*_homozygous") |>
  mutate(has_5x = grepl("5XFAD_carrier", genotype))
plot_pc_grid(pc_plot, "UCI_Trem2-R47H_NSS", "has_5x")

# Counts

counts_trem2 <- make_counts_df(meta_trem2_nss, counts, symbol_map,
                               c("APP", "App", "PSEN1", "Psen1", "Trem2")) |>
  merge(valid_all) |>
  tidyr::pivot_longer(c(APP, App, PSEN1, Psen1, Trem2), names_to = "gene", values_to = "expr")

ggplot(counts_trem2, aes(x = genotype, y = expr, color = valid_expression, shape = genotype)) +
  geom_jitter() +
  facet_wrap(~gene, scales = "free")


# Combine all validation results -----------------------------------------------

valid_samples <- lapply(valid_samples_list, select,
                        any_of(c("unique_specimenID", "valid_variant", "valid_expression"))) |>
  list_rbind() |>
  merge(select(sex_matches, unique_specimenID, valid_sex)) |>
  distinct()

stopifnot(nrow(valid_samples) == nrow(metadata_all))
stopifnot(all(valid_samples$unique_specimenID %in% metadata_all$unique_specimenID))
stopifnot(all(metadata_all$unique_specimenID %in% valid_samples$unique_specimenID))

valid_samples$validated <- valid_samples |>
  select(valid_variant, valid_expression, valid_sex) |>
  mutate(across(everything(), ~is.na(.x) | .x)) |> # NA values count as TRUE
  rowSums() == 3

cat(paste(sum(!valid_samples$validated), "samples failed validation:"), "\n")
print(subset(valid_samples, !validated))
cat("", "\n\n")

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
