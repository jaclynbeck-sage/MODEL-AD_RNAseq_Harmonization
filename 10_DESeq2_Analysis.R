library(synapser)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(stringr)
library(vsn)

source("util_functions.R")

# TODO check diff expression against the genes found in benchmarking that changed
# between reference genomes

# Set up -----------------------------------------------------------------------

file_syn_ids <- config::get("file_syn_ids")
folder_syn_ids <- config::get("staging_syn_ids")
studies <- config::get("studies")

synLogin(silent = TRUE)

github <- paste0(config::get("github_repo_url"),
                 "/blob/main/10_DESeq2_Analysis.R")
tmp_dir <- file.path("output", "tmp")


# Load counts and metadata -----------------------------------------------------

meta_list <- get_all_metadata()
meta_provenance <- do.call(rbind, lapply(meta_list, "[[", "provenance"))

meta_list <- meta_list[studies] |>
  lapply("[[", "data")

symbol_map_file <- synGet(file_syn_ids$symbol_map, downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")
symbol_map <- read.csv(symbol_map_file$path) |>
  arrange(ensembl_gene_id)

counts_list <- get_all_counts_files(studies,
                                    meta_list,
                                    symbol_map,
                                    count_type = "gene_counts")

counts <- do.call(cbind, lapply(counts_list, "[[", "counts")) |>
  # RSEM can return non-integer numbers for counts, so it needs to be rounded
  # to whole numbers for DESeq2
  round(digits = 0)

counts_provenance <- do.call(rbind, lapply(counts_list, "[[", "provenance"))

metadata_all <- do.call(rbind, meta_list) |>
  bin_jax_ages() |>
  mutate(
    sex = str_to_title(sex), # Upper-case
    # Change to "Females" and "Males", plural
    sex_group = paste0(sex, "s"),
    # Add "months" to the end of each age
    age_group = paste(round(ageDeath), "months"),
    # Title case tissue
    tissue = str_to_title(tissue)
  )

# Subset to validated samples only
# TODO temporary until this is finalized and on Synapse
valid_samples <- read.csv("output/Model_AD_valid_samples.csv")
valid_samples <- subset(valid_samples, validated == TRUE)
metadata_all <- subset(metadata_all,
                       unique_specimenID %in% valid_samples$unique_specimenID)

# TODO provenance for valid samples file when it's on Synapse

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, unique_specimenID %in% colnames(counts))
counts <- counts[, metadata_all$unique_specimenID]

provenance_df <- rbind(meta_provenance, counts_provenance) |>
  subset(study %in% studies) |>
  mutate(synID = paste0(id, ".", versionNumber)) |>
  group_by(study) |>
  summarize(used = list(c(synID, file_syn_ids$symbol_map)),
            executed = github) |>
  tibble::column_to_rownames("study")


# Utility functions ------------------------------------------------------------

# Plots the expression of human genes to make sure samples are expressing or
# not expressing them where applicable
plot_human_genes <- function(counts_mat, meta) {
  human_df <- counts_mat[grepl("ENSG", rownames(counts_mat)), ] |>
    as.data.frame() |>
    # Ensembl IDs disappear if we don't make it a column
    tibble::rownames_to_column("ensembl_gene_id") |>
    # All columns except ensembl_gene_id get melted into two columns, sample and
    # count
    tidyr::pivot_longer(cols = where(is.numeric), names_to = "sample",
                        values_to = "count") |>
    # Get genotype information
    merge(meta, by.x = "sample", by.y = "unique_specimenID") |>
    # Get gene symbols
    merge(symbol_map)

  plt <- ggplot(human_df, aes(x = genotype, y = count, color = genotype)) +
    geom_boxplot(position = "dodge", outliers = FALSE) +
    geom_jitter(height = 0) +
    theme_bw() +
    facet_wrap(~gene_symbol, scales = "free", nrow = 1) +
    scale_color_brewer(palette = "Set1")

  print(plt)
}


# Run DESeq2 on a subset of the data. group_vals can either be a string or a
# vector of strings.
run_diff_expr <- function(counts, metadata, group_vals, model = "~ genotype") {
  meta_group <- subset(metadata, group %in% group_vals)
  counts_group <- counts[, meta_group$unique_specimenID]

  # Ensure that any numerical variables that might be in the formula are scaled
  # relative to this group only
  for (col_name in colnames(meta_group)) {
    if (is.numeric(meta_group[, col_name])) {
      meta_group[, col_name] <- scale(meta_group[, col_name])
    }
  }

  # Remove rows with mostly 0 counts -- enough samples have to have a count of
  # 10 or more for each gene
  group_sizes <- meta_group |>
    group_by(genotype, group) |>
    summarize(count = n(), .groups = "drop")

  rows_keep <- rowSums(counts_group >= 10) >= min(group_sizes$count)
  counts_group <- counts_group[rows_keep, ]

  print(str_glue("Subsetting to {sum(rows_keep)} genes."))

  dds <- DESeqDataSetFromMatrix(counts_group, colData = meta_group,
                                design = as.formula(model))

  return(DESeq(dds, quiet = TRUE))
}


verify_sample_sizes <- function(metadata, contrasts, group_val) {
  n_samples <- subset(metadata, group == group_val) |>
    group_by(genotype) |>
    summarize(count = n())

  if (any(n_samples$count < 3)) {
    low_genos <- n_samples$genotype[which(n_samples$count < 3)]
    to_remove <- sapply(contrasts, function(cont) {
      return(cont[2] %in% low_genos | cont[3] %in% low_genos)
    })

    if (any(to_remove)) {
      print(str_glue("----- WARNING for group '{group_val}' -----"))
      print(str_glue("\tThe following genotype(s) do not have enough samples for analysis:"))
      for (item in as.character(low_genos)) {
        print(str_glue("\t\t{item}"))
      }

      print(str_glue("\tThe following contrast(s) will be ignored:"))
      sapply(contrasts[to_remove], function(cont) {
        print(str_glue("\t\t'{cont[2]}' vs '{cont[3]}'"))
      })
      print(str_glue("\n"))
    }
  } else {
    # Nothing should be removed if there are enough samples
    to_remove <- rep(FALSE, length(contrasts))
  }

  return(!to_remove) # vector of contrasts to keep
}


create_model <- function(metadata, group_val, model_vars) {
  meta_group <- subset(metadata, group == group_val)

  to_keep <- sapply(model_vars, function(m_var) {
    if (!(m_var %in% colnames(meta_group))) {
      print(str_glue("----- WARNING for group '{group_val}' -----"))
      print(str_glue("\tVariable '{m_var}' doesn't exist in the metadata and",
                     " will be removed from the model."))
      print(str_glue("\n"))
      return(FALSE)
    }
    if (length(unique(meta_group[, m_var])) < 2) {
      print(str_glue("----- WARNING for group '{group_val}' -----"))
      print(str_glue("\tVariable '{m_var}' has less than 2 unique values. ",
                     "It will be removed from the model."))
      print(str_glue("\n"))
      return(FALSE)
    }
    return(TRUE)
  })

  model <- paste("~", paste(model_vars[to_keep], collapse = " + "))
  return(model)
}


get_all_de_results <- function(metadata, counts, parameters,
                               group_cols = c("sex", "age_group"),
                               model_vars = c("genotype")) {
  meta_sub <- subset(metadata, study == parameters$study) |>
    mutate(genotype = relevel(factor(genotype), ref = parameters$ref_genotype))

  rownames(meta_sub) <- meta_sub$unique_specimenID

  if (length(group_cols) > 1) {
    meta_sub$group <- do.call(paste, meta_sub[, group_cols]) |> factor()
  } else {
    meta_sub$group <- meta_sub[, group_cols] |> factor()
  }

  counts_sub <- counts[, meta_sub$unique_specimenID]

  # Study design
  print(table(meta_sub[, c("genotype", group_cols)]))

  # Human gene expression
  plot_human_genes(counts_sub, meta_sub)

  # Run differential expression on each group separately
  res_all <- lapply(as.character(unique(meta_sub$group)), function(group) {
    to_keep <- verify_sample_sizes(meta_sub, parameters$contrasts, group)
    contrasts <- parameters$contrasts[to_keep]

    if (length(contrasts) == 0) {
      print(str_glue("No contrasts left to compare due to low sample counts. ",
                     "No analysis will be performed for group '{group}'."))
      return(NULL)
    }

    group_model <- create_model(meta_sub, group, model_vars)
    print(str_glue("{group}: using model '{group_model}'"))

    dds <- run_diff_expr(counts_sub, meta_sub,
                         group_vals = group,
                         model = group_model)
    plotDispEsts(dds)

    vsd <- vst(dds, blind = FALSE)
    meanSdPlot(assay(vsd))

    print(plotPCA(vsd, intgroup = "genotype"))

    # Extract results for each contrast
    res_group <- lapply(contrasts, function(contr) {
      res <- results(dds, alpha = 0.05, contrast = contr)
      res <- lfcShrink(dds, res = res, contrast = contr, type = "ashr")

      print(str_glue("Results for group '{group}',"))
      print(str_glue("\t'{contr[2]}' vs '{contr[3]}':"))

      summary(res)

      meta_group <- meta_sub[meta_sub$group == group, ] |>
        select(age_group, sex_group, tissue) |>
        distinct()

      if (nrow(meta_group) > 1) {
        meta_group <- meta_group |>
          summarize(
            across(everything(), ~ paste(sort(unique(.x)), collapse = " & "))
          )
      }

      res <- res |>
        as.data.frame() |>
        tibble::rownames_to_column("ensembl_gene_id") |>
        mutate(model = parameters$model_name,
               case = contr[2],
               control = contr[3],
               age = as.character(meta_group$age_group),
               sex = as.character(meta_group$sex_group),
               tissue = meta_group$tissue) |>
        dplyr::relocate(ensembl_gene_id, .before = baseMean) |>
        dplyr::select(ensembl_gene_id, log2FoldChange, padj, model, case,
                      control, age, sex, tissue)
      return(res)
    })

    return(do.call(rbind, res_group))
  })

  return(do.call(rbind, res_all))
}


get_norm_counts <- function(meta, counts, model_name) {
  counts <- counts[, meta$unique_specimenID]
  sfs <- estimateSizeFactorsForMatrix(counts)

  # Only keep genes that are expressed in at least 3 samples.
  keep <- rowSums(counts > 0) >= 3
  print(str_glue("{model_name}: Keeping {sum(keep)} genes."))

  counts[keep, ] |>
    sageRNAUtils::simple_log2norm(size_factors = sfs) |>
    as.data.frame() |>
    tibble::rownames_to_column("ensembl_gene_id") |>
    tidyr::pivot_longer(cols = -ensembl_gene_id,
                        names_to = "specimenID",
                        values_to = "expression") |>
    merge(select(meta, individualID, unique_specimenID, tissue, sex,
                 age_group, genotype),
          by.x = "specimenID", by.y = "unique_specimenID") |>
    mutate(age = age_group,
           model = model_name) |>
    # Put columns in a specific order
    dplyr::select(ensembl_gene_id, individualID, expression, tissue, sex, age,
                  genotype, model)
}


# Differential expression for individual studies -------------------------------

## Jax.IU.Pitt_5XFAD -----------------------------------------------------------

# This data has no separate batches.
meta_jax5x <- subset(metadata_all, study == "Jax.IU.Pitt_5XFAD") |>
  subset(age_group != "8 months") # We only want 4 mo and 12 mo for the explorer

params_jax5x <- list(
  study = "Jax.IU.Pitt_5XFAD",
  model_name = "5xFAD (IU/Jax/Pitt)",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier"))
)

# Separated by sex and ageDeath
res_jax5x <- get_all_de_results(
  meta_jax5x, counts, params_jax5x,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype")
)

# Males and females together, separated by ageDeath
res_jax5x_mf <- get_all_de_results(
  meta_jax5x, counts, params_jax5x,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex")
)

jax5x_de_file <- str_glue("output/de_output/{params_jax5x$study}_differential_expression.csv")
write.csv(rbind(res_jax5x, res_jax5x_mf), jax5x_de_file,
          row.names = FALSE, quote = FALSE)

norm_jax5x <- get_norm_counts(meta_jax5x, counts, params_jax5x$model_name)
jax5x_norm_file <- str_glue("output/de_output/{params_jax5x$study}_normalized_expression.csv")
write.csv(norm_jax5x, jax5x_norm_file,
          row.names = FALSE, quote = FALSE)

# Upload to Synapse
synapse_upload_de(params_jax5x$study, provenance_df, folder_syn_ids,
                  jax5x_de_file, jax5x_norm_file)


## Jax.IU.Pitt_APOE4.Trem2.R47H ------------------------------------------------

meta_load1 <- subset(metadata_all, study == "Jax.IU.Pitt_APOE4.Trem2.R47H") |>
  # drop the two samples with Trem2-R47H_heterozygous genotypes
  subset(!grepl("heterozygous", genotype)) |>
  # We only want 4, 12, and 24 months for the explorer
  subset(age_group != "8 months")

ref_geno <- "APOE4-KI_WT; Trem2-R47H_WT"

params_load1 <- list(
  study = "Jax.IU.Pitt_APOE4.Trem2.R47H",
  model_name = "LOAD1",
  ref_genotype = ref_geno,
  contrasts = list(
    c("genotype", "APOE4-KI_homozygous; Trem2-R47H_homozygous", ref_geno),
    c("genotype", "APOE4-KI_homozygous; Trem2-R47H_WT", ref_geno),
    c("genotype", "APOE4-KI_WT; Trem2-R47H_homozygous", ref_geno)
  )
)

# The ages of the mice span a wide range within each bin. The 4 month age
# group spans ages 3-4 months, and the 12 month group actually spans ages
# 13-16 months. Because of this spread, we will use their real age as a
# variable in the model within each age group.

# Separated by sex and ageDeath
res_load1 <- get_all_de_results(
  meta_load1, counts, params_load1,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype", "ageDeathNumeric", "sequencingBatch")
)

# Males and females together, separated by ageDeath
res_load1_mf <- get_all_de_results(
  meta_load1, counts, params_load1,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex", "ageDeathNumeric", "sequencingBatch")
)

res_load1_all <- rbind(res_load1, res_load1_mf) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(case,
                            "APOE4-KI_WT; Trem2-R47H_homozygous" ~ "Trem2R47H",
                            "APOE4-KI_homozygous; Trem2-R47H_WT" ~ "APOE4",
                            .default = model))

load1_de_file <- str_glue("output/de_output/{params_load1$study}_differential_expression.csv")
write.csv(res_load1_all, load1_de_file,
          row.names = FALSE, quote = FALSE)

norm_load1 <- get_norm_counts(meta_load1, counts, params_load1$model_name) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(genotype,
                            "APOE4-KI_WT; Trem2-R47H_homozygous" ~ "Trem2R47H",
                            "APOE4-KI_homozygous; Trem2-R47H_WT" ~ "APOE4",
                            .default = model))

# Split into one file per "model" for the explorer
load1_norm_files <- lapply(unique(norm_load1$model), function(model_name) {
  norm_data <- subset(norm_load1, model == model_name | genotype == ref_geno) |>
    mutate(model = model_name) # Fix WT model to match

  output_name <- case_match(model_name,
                            "LOAD1" ~ "APOE4.Trem2.R47H",
                            "Trem2R47H" ~ "Trem2.R47H",
                            "APOE4" ~ "APOE4")

  norm_filename <- str_glue("output/de_output/Jax.IU.Pitt_{output_name}_normalized_expression.csv")
  write.csv(norm_data, norm_filename,
            row.names = FALSE, quote = FALSE)

  return(norm_filename)
})

# Upload to Synapse
synapse_upload_de(params_load1$study, provenance_df, folder_syn_ids,
                  load1_de_file, load1_norm_files)


## Jax.IU.Pitt_LOAD2 -----------------------------------------------------------

meta_load2 <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2") |>
  # Drop the one hAPP-3/3 heterozygous genotype.
  subset(!grepl("heterozygous", genotype)) |>
  # Note: for the explorer, we only want to compare LOAD1 vs LOAD2. There is
  # only one age group (18 months) with enough WT samples to compare LOAD2 vs
  # WT, so this analysis is skipped. We also remove WT mice from the data set,
  # because either:
  #   - No WT mice exist in some age/sex groups already, or
  #   - In the non-18-month age groups where WT mice exist, LOAD1/LOAD2 were
  #     sequenced in two batches together but all WT mice were sequenced
  #     separately in a third batch, which doesn't allow us to use
  #     sequencingBatch in the model.
  subset(genotype != "APOE4-KI_WT; Trem2-R47H_WT; hAPP-3/3_WT")

load2_geno <- "APOE4-KI_homozygous; Trem2-R47H_homozygous; hAPP-3/3_homozygous"
load1_geno <- "APOE4-KI_homozygous; Trem2-R47H_homozygous; hAPP-3/3_WT"

params_load2 <- list(
  study = "Jax.IU.Pitt_LOAD2",
  model_name = "LOAD2",
  ref_genotype = load1_geno,
  contrasts = list(c("genotype", load2_geno, load1_geno))
)

# The ages of the mice span a wide range within each bin. The 4 month age
# group spans ages 3-6 months, and the 12 month group spans ages 12-13 months,
# the 18 month group spans ages 18-19 months, and the 24 month group spans
# ages 22-25 months. Because of this spread, we will use their real age as a
# variable in the model within each age group. Some age/sex groups have two
# batches, so we also use sequencingBatch in the model for these groups.

# Separated by sex and ageDeath
res_load2 <- get_all_de_results(
  meta_load2, counts, params_load2,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype", "ageDeathNumeric", "sequencingBatch")
)

# Males and females together, separated by ageDeath
res_load2_mf <- get_all_de_results(
  meta_load2, counts, params_load2,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex", "ageDeathNumeric", "sequencingBatch")
)

load2_de_file <- str_glue("output/de_output/{params_load2$study}_differential_expression.csv")
write.csv(rbind(res_load2, res_load2_mf), load2_de_file,
          row.names = FALSE, quote = FALSE)

norm_load2 <- get_norm_counts(meta_load2, counts, params_load2$model_name)
load2_norm_file <- str_glue("output/de_output/{params_load2$study}_normalized_expression.csv")
write.csv(norm_load2, load2_norm_file,
          row.names = FALSE, quote = FALSE)


# Upload to Synapse
synapse_upload_de(params_load2$study, provenance_df, folder_syn_ids,
                  load2_de_file, load2_norm_file)


## Jax.IU.Pitt_LOAD2.PrimaryScreen ---------------------------------------------

# This study is actually multiple different models together in one file, so we
# break the models out here. This study has 7 sequencing batches with the
# following splits:
#   21-model-ad-007: LOAD2.Bin1-KI_KI/KI and LOAD2.Bin1-KI_WT/WT, and
#                    LOAD2.Ptprb-D57N_KI/KI and LOAD2.Ptprb-D57N_WT/WT
#   22-model-ad-003: LOAD2.Epha1-KI_KI/KI and LOAD2.Epha1-KI_WT/WT, some C57BL6J
#   22-model-ad-007: LOAD2.Cd2ap_KI/KI and LOAD2.Cd2ap_WT/WT, and
#                    LOAD2.hTau-H2_homozygous and LOAD2.hTau-H2_WT, some C57BL6J
#   23-model-ad-001: some LOAD2.IL1rap-Exon3KO_homozygous and LOAD2.IL1rap-Exon3KO_WT, some C57BL6J
#   23-model-ad-003: some LOAD2, some LOAD2.Adamts4-KO_homozygous,
#                    some LOAD2.Il34-Y213_homozygous, some LOAD2.Ptk2b-intron5SNP_homozygous,
#                    some LOAD2.Scimp-upstreamSNP_homozygous
#   23-model-ad-003-run2: some LOAD2, some LOAD2.Adamts4-KO_homozygous and LOAD2.Adamts4-KO_WT,
#                    some LOAD2.Il34-Y213_homozygous, some LOAD2.Ptk2b-intron5SNP_homozygous,
#                    some LOAD2.Scimp-upstreamSNP_homozygous
#   24-model-ad-003: some LOAD2.IL1rap-Exon3KO_homozygous and LOAD2.IL1rap-Exon3KO_WT
#
# Genotype pairings:
#   LOAD2.Bin1-KI_KI/KI and LOAD2.Bin1-KI_WT/WT
#     - all in the same batch
#     - 5-6 KI and 3 WT per group
#   LOAD2.Ptprb-D57N_KI/KI and LOAD2.Ptprb-D57N_WT/WT
#     - all in the same batch
#     - 6 KI and 3 WT per group
#   LOAD2.Epha1-KI_KI/KI and LOAD2.Epha1-KI_WT/WT
#     - all in the same batch
#     - 6 KI and 6 WT per group
#   LOAD2.Cd2ap_KI/KI and LOAD2.Cd2ap_WT/WT
#     - all in the same batch
#     - 6 KI and 3 WT per group
#   LOAD2.hTau-H2_homozygous and LOAD2.hTau-H2_WT
#     - all in the same batch
#     - 6 homozygous and 2-3 WT per group
#       * 12 month male WT only has 2 samples and will be excluded from analysis.
#   LOAD2.IL1rap-Exon3KO_homozygous and LOAD2.IL1rap-Exon3KO_WT
#     - split into two batches:
#       - all 4 month samples are in one batch, and 12 month samples are split between the two
#       - sexes and genotypes are split evenly across batches
#       * 12 month samples will need batch correction
#     - 4 months: 5-6 homozygous and 3 WT per group
#     - 12 months: 26 homozygous and 21-23 WT per group
#       - differences in diet and coming from two centers. This info needs to be added to the metadata file
#   LOAD2.Adamts4-KO_homozygous and LOAD2.Adamts4-KO_WT (and LOAD2?)
#     - split into two batches:
#       * 23-model-ad-003 should be thrown out. There is only 1 sample for each group + genotype in it
#     - LOAD2.Adamts4-KO_WT only has 2 samples total. Unclear why they are labeled separately from LOAD2
#       since these should be equivalent
#     - 4 months: 5 homozygous and 5-6 LOAD2 per group
#     - 12 month males: 5 homozygous and 6 LOAD2
#     - 12 month females: 6 homozygous, 2 LOAD2, and 2 LOAD2.Adamts4-KO_WT
#   LOAD2.Il34-Y213_homozygous and LOAD2
#     - split into 2 batches, 23-model-ad-003 should be thrown out
#     - 4 months: 5-6 homozygous and 5-6 LOAD2 per group
#     - 12 month males: 5 homozygous and 6 LOAD2
#     - 12 month females: 6 homozygous and 2 LOAD2
#       * This group should get thrown out
#   LOAD2.Ptk2b-intron5SNP_homozygous and LOAD2
#     - split into 2 batches, 23-model-ad-003 should be thrown out
#     - 4 months: 5-6 homozygous and 5-6 LOAD2 per group
#     - 12 month males: 4 homozygous and 6 LOAD2
#     - 12 month females: 4 homozygous and 2 LOAD2
#       * This group should get thrown out
#   LOAD2.Scimp-upstreamSNP_homozygous and LOAD2
#     - split into 2 batches, 23-model-ad-003 should be thrown out
#     - 4 months: 5-6 homozygous and 5-6 LOAD2 per group
#     - 12 month males: 6 homozygous and 6 LOAD2
#     - 12 month females: 6 homozygous and 2 LOAD2
#       * This group should get thrown out
#
# Note: Adamts4, Il34, Ptk2b, and Scimp were all bred at IU, which is why their
# batches/designs are different
#
# Note: real ages vary by up to 3 months, so the numeric age needs to be a
# covariate in the design for DEseq2.

meta_load2pri <- subset(metadata_all, study == "Jax.IU.Pitt_LOAD2.PrimaryScreen") |>
  mutate(
    model = str_replace(genotype, "_(homozygous|WT|KI/KI)(/WT)?", "")
  )

# TODO


## UCI_3xTg-AD -----------------------------------------------------------------

params_3x <- list(
  study = "UCI_3xTg-AD",
  model_name = "3xTg-AD",
  ref_genotype = "3xTg-AD_noncarrier",
  contrasts = list(c("genotype", "3xTg-AD_carrier", "3xTg-AD_noncarrier"))
)

# There are no batch effects, no separate tissues
res_3x <- get_all_de_results(
  metadata_all, counts, params_3x,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype")
)

res_3x_mf <- get_all_de_results(
  metadata_all, counts, params_3x,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex")
)

de_3x_file <- str_glue("output/de_output/{params_3x$study}_differential_expression.csv")
write.csv(rbind(res_3x, res_3x_mf), de_3x_file,
          row.names = FALSE, quote = FALSE)

norm_3x <- get_norm_counts(subset(metadata_all, study == params_3x$study),
                           counts, params_3x$model_name)

norm_3x_file <- str_glue("output/de_output/{params_3x$study}_normalized_expression.csv")
write.csv(norm_3x, norm_3x_file,
          row.names = FALSE, quote = FALSE)

# Upload to Synapse
synapse_upload_de(params_3x$study, provenance_df, folder_syn_ids,
                  de_3x_file, norm_3x_file)


## UCI_5XFAD -------------------------------------------------------------------

# This data has a sequencing batch, but it corresponds to the age of death so
# we don't include it in the model. There are two tissues.

# Drop 8 month time point
meta_uci5x <- subset(metadata_all,
                     study == "UCI_5XFAD" & age_group != "8 months")

params_5x <- list(
  study = "UCI_5XFAD",
  model_name = "5xFAD (UCI)",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier"))
)

res_5x <- get_all_de_results(
  meta_uci5x, counts, params_5x,
  group_cols = c("tissue", "sex", "age_group"),
  model_vars = c("genotype")
)

res_5x_mf <- get_all_de_results(
  meta_uci5x, counts, params_5x,
  group_cols = c("tissue", "age_group"),
  model_vars = c("genotype", "sex")
)

de_5x_file <- str_glue("output/de_output/{params_5x$study}_differential_expression.csv")
write.csv(rbind(res_5x, res_5x_mf), de_5x_file,
          row.names = FALSE, quote = FALSE)

norm_5x <- get_norm_counts(meta_uci5x, counts, params_5x$model_name)
norm_5x_file <- str_glue("output/de_output/{params_5x$study}_normalized_expression.csv")
write.csv(norm_5x, norm_5x_file,
          row.names = FALSE, quote = FALSE)

# Upload to Synapse
synapse_upload_de(params_5x$study, provenance_df, folder_syn_ids,
                  de_5x_file, norm_5x_file)


## UCI_ABCA7 -------------------------------------------------------------------

# This data has rnaBatch with a bad split:
#   'Abca7-V1599M_homozygous' and '5XFAD_carrier; Abca7-V1599M_homozygous' were
#     extracted together across two batches
#   `5XFAD_carrier` and `5XFAD_noncarrier` were extracted together across two
#     other batches
#   Our comparisons are between genotypes with no batch overlap
# rnaBatch also has exact overlap with age of death, so each age of death contains
# one batch per pair of genotypes.
# Given this, we can't actually batch correct this data, which is not ideal.

params_abca7 <- list(
  study = "UCI_ABCA7",
  model_name = "Abca7*V1599M",
  ref_genotype = "5XFAD_noncarrier",
  # We want Abca7-homozygous vs WT, and Abca7-5xFAD vs 5xFAD for the explorer
  contrasts = list(c("genotype", "Abca7-V1599M_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Abca7-V1599M_homozygous", "5XFAD_carrier"))
)

res_abca7 <- get_all_de_results(
  metadata_all, counts, params_abca7,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype")
)

res_abca7_mf <- get_all_de_results(
  metadata_all, counts, params_abca7,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex")
)

res_abca7_all <- rbind(res_abca7, res_abca7_mf) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = ifelse(case == "5XFAD_carrier; Abca7-V1599M_homozygous",
                        "Abca7*V1599M.5xFAD", model))

de_abca7_file <- str_glue("output/de_output/{params_abca7$study}_differential_expression.csv")
write.csv(res_abca7_all, de_abca7_file,
          row.names = FALSE, quote = FALSE)

norm_abca7 <- get_norm_counts(subset(metadata_all, study == params_abca7$study),
                              counts, params_abca7$model_name) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(genotype,
                            "5XFAD_carrier; Abca7-V1599M_homozygous" ~ "Abca7*V1599M.5xFAD",
                            "5XFAD_carrier" ~ "Abca7*V1599M.5xFAD",
                            .default = model))

# Split into one file per "model" for the explorer
norm_abca7_files <- lapply(unique(norm_abca7$model), function(model_name) {
  norm_data <- subset(norm_abca7, model == model_name)

  output_name <- str_replace(model_name, "\\*V1599M", "") |> str_to_upper()

  norm_file <- str_glue("output/de_output/UCI_{output_name}_normalized_expression.csv")
  write.csv(norm_data, norm_file,
            row.names = FALSE, quote = FALSE)
  return(norm_file)
})

# Upload to Synapse
synapse_upload_de(params_abca7$study, provenance_df, folder_syn_ids,
                  de_abca7_file, norm_abca7_files)


## UCI_Bin1K358R ---------------------------------------------------------------

# This study has some non-ideal batch splits. For each age group, 5XFAD_carrier
# and 5XFAD_noncarrier are in one rnaBatch, while Bin1-K358R_homozygous and
# 5XFAD_carrier; Bin1-K358R_homozygous are in another. These are not the pairs
# of comparisons we want to make for DE, so we can't batch correct for this
# variable in a fixed-effect model. A mixed-effect model may be more appropriate.
#
# Each age group has two `libraryBatch`es, with each genotype split fairly
# evenly between the two batches. However, means that about half of the batches
# have < 3 samples per genotype in them. There is also a typo batch with exactly
# 1 sample (06.0.22) which was probably intended to be "06.01.22". The split
# between Dec-22 and Nov-22 is suspiciously uneven, with either 1 sample or
# (N-1) samples in each batch for each genotype. Given all of this, we can not
# batch correct for this variable.
#
# `sequencingBatch` corresponds exactly to age group and can be ignored.
#
# Genotype differences (5XFAD_carrier vs 5XFAD_noncarrier) appear to be driving
# the first two PCAs far more than rnaBatch, however some separation by rnaBatch
# is visible for some age groups.

# PCAs to look at batch
meta_bin1 <- subset(metadata_all, study == "UCI_Bin1K358R")
counts_bin1 <- counts[, meta_bin1$unique_specimenID]
cpm_bin1 <- sweep(counts_bin1, 2, colSums(counts_bin1), "/") * 1e6
log_bin1 <- log2(cpm_bin1 + 0.5)

do_plots <- function(pc_df) {
  plt <- ggplot(pc_df, aes(x = PC1, y = PC2, color = genotype, shape = rnaBatch)) +
    geom_point() + theme_bw() +
    facet_wrap(~age_group)
  print(plt)

  plt <- ggplot(pc_df, aes(x = PC1, y = PC2, color = rnaBatch, shape = genotype)) +
    geom_point() + theme_bw() +
    facet_wrap(~age_group)
  print(plt)

  plt <- ggplot(pc_df, aes(x = PC1, y = PC2, color = genotype)) +
    geom_point() + theme_bw() +
    facet_grid(rows = vars(age_group), cols = vars(rnaBatch))
  print(plt)

  plt <- ggplot(pc_df, aes(x = PC1, y = PC2, color = rnaBatch)) +
    geom_point() + theme_bw() +
    facet_grid(rows = vars(age_group), cols = vars(genotype))
  print(plt)
}

# All samples together -- primary driver is whether there is a 5X or not
pc <- prcomp(t(log_bin1))
pc$x <- merge(pc$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc$x)

# 4 month age group -- primary driver is whether there is a 5X or not
pc_age4 <- prcomp(t(log_bin1[, meta_bin1$age_group == "4 months"]))
pc_age4$x <- merge(pc_age4$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_age4$x)

# 12 month age group -- primary driver is whether there is a 5X or not
pc_age12 <- prcomp(t(log_bin1[, meta_bin1$age_group == "12 months"]))
pc_age12$x <- merge(pc_age12$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_age12$x)

# Bin1 vs WT -- primary driver is age group
pc_b1 <- prcomp(t(log_bin1[, meta_bin1$genotype %in% c("Bin1-K358R_homozygous", "5XFAD_noncarrier")]))
pc_b1$x <- merge(pc_b1$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b1$x)

# 5x.Bin1 vs 5x -- primary driver is age group
pc_b5x <- prcomp(t(log_bin1[, meta_bin1$genotype %in% c("5XFAD_carrier; Bin1-K358R_homozygous", "5XFAD_carrier")]))
pc_b5x$x <- merge(pc_b5x$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b5x$x)

# Bin1 vs WT 4 months -- no visual difference
pc_b4 <- prcomp(t(log_bin1[, meta_bin1$age_group == "4 months" &
                             meta_bin1$genotype %in% c("Bin1-K358R_homozygous", "5XFAD_noncarrier")]))
pc_b4$x <- merge(pc_b4$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b4$x)

# Bin1 vs WT 12 months -- slight separation by genotype / batch
pc_b12 <- prcomp(t(log_bin1[, meta_bin1$age_group == "12 months" &
                             meta_bin1$genotype %in% c("Bin1-K358R_homozygous", "5XFAD_noncarrier")]))
pc_b12$x <- merge(pc_b12$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b12$x)

# 5x.Bin1 vs 5x 4 months -- no visual difference
pc_b5x4 <- prcomp(t(log_bin1[, meta_bin1$age_group == "4 months" &
                              meta_bin1$genotype %in% c("5XFAD_carrier; Bin1-K358R_homozygous", "5XFAD_carrier")]))
pc_b5x4$x <- merge(pc_b5x4$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b5x4$x)

# 5x.Bin1 vs 5x 12 months -- slight separation by genotype / batch
pc_b5x12 <- prcomp(t(log_bin1[, meta_bin1$age_group == "12 months" &
                               meta_bin1$genotype %in% c("5XFAD_carrier; Bin1-K358R_homozygous", "5XFAD_carrier")]))
pc_b5x12$x <- merge(pc_b5x12$x, meta_bin1, by.x = "row.names", by.y = "unique_specimenID")

do_plots(pc_b5x12$x)


params_bin1 <- list(
  study = "UCI_Bin1K358R",
  model_name = "Bin1-K358R",
  ref_genotype = "5XFAD_noncarrier",
  # We want Bin1-K358R_homozygous vs WT, and Bin1-5xFAD vs 5xFAD for the explorer
  contrasts = list(c("genotype", "Bin1-K358R_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Bin1-K358R_homozygous", "5XFAD_carrier"))
)

res_bin1 <- get_all_de_results(
  metadata_all, counts, params_bin1,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype")
)

res_bin1_mf <- get_all_de_results(
  metadata_all, counts, params_bin1,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex")
)

res_bin1_all <- rbind(res_bin1, res_bin1_mf) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = ifelse(case == "5XFAD_carrier; Bin1-K358R_homozygous",
                        "Bin1-K358R.5xFAD", model))

de_bin1_file <- str_glue("output/de_output/{params_bin1$study}_differential_expression.csv")
write.csv(res_bin1_all, de_bin1_file,
          row.names = FALSE, quote = FALSE)

norm_bin1 <- get_norm_counts(subset(metadata_all, study == params_bin1$study),
                             counts, params_bin1$model_name) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(genotype,
                            "5XFAD_carrier; Bin1-K358R_homozygous" ~ "Bin1-K358R.5xFAD",
                            "5XFAD_carrier" ~ "Bin1-K358R.5xFAD",
                            .default = model))

# Split into one file per "model" for the explorer
norm_bin1_files <- lapply(unique(norm_bin1$model), function(model_name) {
  norm_data <- subset(norm_bin1, model == model_name)

  norm_file <- str_glue("output/de_output/UCI_{model_name}_normalized_expression.csv")
  write.csv(norm_data, norm_file,
            row.names = FALSE, quote = FALSE)
  return(norm_file)
})

# Upload to Synapse
synapse_upload_de(params_bin1$study, provenance_df, folder_syn_ids,
                  de_bin1_file, norm_bin1_files)


## UCI_Clu-h2kbKI --------------------------------------------------------------

# For this study, there are two sequencing batches but they correspond to age
# group: all 4 month samples are in one batch and all 12 month samples are in
# another. Therefore we do not need to batch correct the data.

# This study has two tissues.

params_clu <- list(
  study = "UCI_Clu-h2kbKI",
  model_name = "Clu-h2kbKI",
  ref_genotype = "5XFAD_noncarrier",
  # We want Clu-homozygous vs WT, and Clu-5xFAD vs 5xFAD for the explorer
  contrasts = list(c("genotype", "Clu-rs2279590_KI_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Clu-rs2279590_KI_homozygous", "5XFAD_carrier"))
)

res_clu <- get_all_de_results(
  metadata_all, counts, params_clu,
  group_cols = c("tissue", "sex", "age_group"),
  model_vars = c("genotype")
)

res_clu_mf <- get_all_de_results(
  metadata_all, counts, params_clu,
  group_cols = c("tissue", "age_group"),
  model_vars = c("genotype", "sex")
)

res_clu_all <- rbind(res_clu, res_clu_mf) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = ifelse(case == "5XFAD_carrier; Clu-rs2279590_KI_homozygous",
                        "Clu-h2kbKI.5xFAD", model))

de_clu_file <- str_glue("output/de_output/{params_clu$study}_differential_expression.csv")
write.csv(res_clu_all, de_clu_file,
          row.names = FALSE, quote = FALSE)

norm_clu <- get_norm_counts(subset(metadata_all, study == params_clu$study),
                            counts, params_clu$model_name) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(genotype,
                            "5XFAD_carrier; Clu-rs2279590_KI_homozygous" ~ "Clu-h2kbKI.5xFAD",
                            "5XFAD_carrier" ~ "Clu-h2kbKI.5xFAD",
                            .default = model))

# Split into one file per "model" for the explorer
norm_clu_files <- lapply(unique(norm_clu$model), function(model_name) {
  norm_data <- subset(norm_clu, model == model_name)

  norm_file <- str_glue("output/de_output/UCI_{model_name}_normalized_expression.csv")
  write.csv(norm_data, norm_file,
            row.names = FALSE, quote = FALSE)
  return(norm_file)
})

# Upload to Synapse
synapse_upload_de(params_clu$study, provenance_df, folder_syn_ids,
                  de_clu_file, norm_clu_files)


## UCI_Trem2-R47H_NSS ----------------------------------------------------------

# This data has all 3 batch variables filled in (libraryBatch, sequencingBatch,
# rnaBatch). None of them are entirely unique combinations of another even when
# split by age and sex. All 3 batch variables have the same issues that the
# Abca7 study does, where the genotypes we want to compare were extracted /
# prepped / sequenced in separate batches.
# - "libraryBatch" for the 12 month time point doesn't have this issue, however,
#   the batches are only 2-3 mice.
# Given this, we can't actually batch correct this data, which is not ideal.

params_trem2nss <- list(
  study = "UCI_Trem2-R47H_NSS",
  model_name = "Trem2-R47H_NSS",
  ref_genotype = "5XFAD_noncarrier",
  # We want Trem2R47H vs WT, and Trem2R47H-5xFAD vs 5xFAD for the explorer
  contrasts = list(c("genotype", "Trem2-R47H_NSS_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Trem2-R47H_NSS_homozygous", "5XFAD_carrier"))
)

res_trem2nss <- get_all_de_results(
  metadata_all, counts, params_trem2nss,
  group_cols = c("sex", "age_group"),
  model_vars = c("genotype")
)

res_trem2nss_mf <- get_all_de_results(
  metadata_all, counts, params_trem2nss,
  group_cols = c("age_group"),
  model_vars = c("genotype", "sex")
)

res_trem2nss_all <- rbind(res_trem2nss, res_trem2nss_mf) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = ifelse(case == "5XFAD_carrier; Trem2-R47H_NSS_homozygous",
                        "Trem2-R47H_NSS.5xFAD", model))

de_trem2nss_file <- str_glue("output/de_output/{params_trem2nss$study}_differential_expression.csv")
write.csv(res_trem2nss_all, de_trem2nss_file,
          row.names = FALSE, quote = FALSE)

norm_trem2nss <- get_norm_counts(subset(metadata_all, study == params_trem2nss$study),
                                 counts, params_trem2nss$model_name) |>
  # Re-name model for certain genotypes to work with the explorer
  mutate(model = case_match(genotype,
                            "5XFAD_carrier; Trem2-R47H_NSS_homozygous" ~ "Trem2-R47H_NSS.5xFAD",
                            "5XFAD_carrier" ~ "Trem2-R47H_NSS.5xFAD",
                            .default = model))

# Split into one file per "model" for the explorer
norm_trem2nss_files <- lapply(unique(norm_trem2nss$model), function(model_name) {
  norm_data <- subset(norm_trem2nss, model == model_name)

  output_name <- str_replace(model_name, "5xFAD", "5XFAD")

  norm_file <- str_glue("output/de_output/UCI_{output_name}_normalized_expression.csv")
  write.csv(norm_data, norm_file,
            row.names = FALSE, quote = FALSE)

  return(norm_file)
})

# Upload to Synapse
synapse_upload_de(params_trem2nss$study, provenance_df, folder_syn_ids,
                  de_trem2nss_file, norm_trem2nss_files)
