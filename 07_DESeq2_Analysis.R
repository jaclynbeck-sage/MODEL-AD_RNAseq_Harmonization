library(synapser)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(stringr)
library(vsn)

# TODO check diff expression against the genes found in benchmarking that changed
# between reference genomes

# Set up -----------------------------------------------------------------------

# Synapse IDs used in this script -- merged metadata and the gene_counts
# file with data from all studies
syn_metadata_file_id <- "syn61850266"
syn_gene_counts_id <- "syn62690577"
syn_symbol_map_id <- "syn62063692"

synLogin()

github <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/07_DESeq2_Analysis.Rmd"
tmp_dir <- file.path("data", "tmp")


# Load counts and metadata -----------------------------------------------------

metadata_file <- synGet(syn_metadata_file_id, downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")
counts_file <- synGet(syn_gene_counts_id, downloadLocation = tmp_dir,
                      ifcollision = "overwrite.local")
symbol_map_file <- synGet(syn_symbol_map_id, downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

metadata_all <- read.csv(metadata_file$path)

counts <- read.table(counts_file$path, header = TRUE, sep = "\t",
                     row.names = 1) %>%
  # Get rid of transcript_id column
  select(-transcript_id.s.) %>%
  # RSEM can return non-integer numbers for counts, so it needs to be rounded
  # to whole numbers for DESeq2
  round(digits = 0)

# Subset to validated samples only
# TODO temporary until this is finalized and on Synapse
valid_samples <- read.csv("data/Model_AD_valid_samples.csv")
valid_samples <- subset(valid_samples, validated == TRUE)
metadata_all <- subset(metadata_all,
                       merged_file_specimenID %in% valid_samples$merged_file_specimenID)

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, merged_file_specimenID %in% colnames(counts))
counts <- counts[, metadata_all$merged_file_specimenID]

symbol_map <- read.csv(symbol_map_file$path)


# Utility functions ------------------------------------------------------------

# Plots the expression of human genes to make sure samples are expressing or
# not expressing them where applicable
plot_human_genes <- function(counts_mat, meta) {
  human_df <- counts_mat[grepl("ENSG", rownames(counts_mat)), ] %>%
    # Ensembl IDs disappear if we don't make it a column
    mutate(ensembl_gene_id = rownames(.)) %>%
    # All columns except ensembl_gene_id get melted into two columns, sample and
    # count
    tidyr::pivot_longer(cols = where(is.numeric), names_to = "sample",
                        values_to = "count") %>%
    # Get genotype information
    merge(meta, by.x = "sample", by.y = "merged_file_specimenID") %>%
    # Get gene symbols
    merge(symbol_map)

  plt <- ggplot(human_df, aes(x = genotype, y = count, color = genotype)) +
    geom_boxplot(position = "dodge", outliers = FALSE) +
    geom_jitter(height = 0) +
    theme_bw() +
    facet_wrap(~gene_symbol, scales = "free") +
    scale_color_brewer(palette = "Set1")

  print(plt)
}


# Run DESeq2 on a subset of the data. group_vals can either be a string or a
# vector of strings.
run_diff_expr <- function(counts, metadata, group_vals, model = "~ genotype") {
  meta_group <- subset(metadata, group %in% group_vals)
  counts_group <- counts[, meta_group$merged_file_specimenID]

  # Ensure that any numerical variables that might be in the formula are scaled
  # relative to this group only
  for (col_name in colnames(meta_group)) {
    if (is.numeric(meta_group[, col_name])) {
      meta_group[, col_name] <- scale(meta_group[, col_name])
    }
  }

  # Remove rows with mostly 0 counts -- enough samples have to have a count of
  # 10 or more for each gene
  group_sizes <- meta_group %>%
    group_by(genotype, group) %>%
    summarize(count = n(), .groups = "drop")

  rows_keep <- rowSums(counts_group >= 10) >= min(group_sizes$count)
  counts_group <- counts_group[rows_keep, ]

  print(str_glue("Subsetting to {sum(rows_keep)} genes."))

  dds <- DESeqDataSetFromMatrix(counts_group, colData = meta_group,
                                design = as.formula(model))

  return(DESeq(dds, quiet = TRUE))
}


verify_sample_sizes <- function(metadata, contrasts, group_val) {
  n_samples <- subset(metadata, group == group_val) %>%
    group_by(genotype) %>%
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


get_all_de_results <- function(metadata_all, counts, study, ref_genotype,
                               contrasts, group_cols = c("sex", "ageDeath"),
                               model_vars = c("genotype")) {
  meta_sub <- subset(metadata_all, study_name == study) %>%
    mutate(genotype = relevel(factor(genotype), ref = ref_genotype))

  meta_sub$group <- do.call(paste, meta_sub[, group_cols]) %>%
    paste("months") %>% factor()

  counts_sub <- counts[, meta_sub$merged_file_specimenID]

  # Study design
  print(table(meta_sub[, c("genotype", group_cols)]))

  # Human gene expression
  plot_human_genes(counts_sub, meta_sub)

  # Run differential expression on each group separately
  res_all <- lapply(unique(meta_sub$group), function(group) {
    to_keep <- verify_sample_sizes(meta_sub, contrasts, group)
    contrasts <- contrasts[to_keep]

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

      res <- res %>%
        as.data.frame() %>%
        mutate(ensembl_gene_id = rownames(.),
               group = group,
               test_genotype = contr[2],
               reference_genotype = contr[3],
               contrast = str_glue("{contr[2]} - {contr[3]}"))
      return(res)
    })

    return(list(dds_result = do.call(rbind, res_group),
                norm_counts = assay(vsd)))
  })

  norm_counts <- lapply(res_all, "[[", "norm_counts")
  names(norm_counts) <- make.names(paste(unique(meta_sub$group)))

  return(list(dds_result = do.call(rbind, lapply(res_all, "[[", "dds_result")),
              norm_counts = res_all$norm_counts))
}


# Differential expression for individual studies -------------------------------

## Jax.IU.Pitt_5XFAD -----------------------------------------------------------

# This data has no separate batches.
# Bin the ages into 4, 6, and 12 months
meta_jax5x <- subset(metadata_all, study_name == "Jax.IU.Pitt_5XFAD") %>%
  mutate(ageDeath = case_when(ageDeath < 5 ~ 4,
                              ageDeath > 5 & ageDeath < 10 ~ 6,
                              ageDeath > 10 ~ 12))

res_jax5x <- get_all_de_results(
  meta_jax5x, counts, study = "Jax.IU.Pitt_5XFAD",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier")),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("genotype")
)
# TODO write to file


## Jax.IU.Pitt_APOE4.Trem2.R47H ------------------------------------------------

meta_load1 <- subset(metadata_all, study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H") %>%
  # drop the two samples with Trem2-R47H_heterozygous genotypes, drop any
  # samples with NA genotypes
  subset(!grepl("heterozygous", genotype) & !is.na(genotype)) %>%

  # The ages of the mice span a wide range and need to be binned into 4, 8, 12,
  # and 24 month groups. The 4 month age group spans 3-4 months, the 8 month
  # group spans 7-10 months, the 12 month group actually spans 13-16 months, and
  # the 24 month group spans 24-28 months. Because of this spread, we will use
  # their real age as a variable in the model within each age group.
  mutate(real_age = ageDeath,
         ageDeath = case_when(ageDeath < 5 ~ 4,
                              ageDeath > 5 & ageDeath < 11 ~ 8,
                              ageDeath > 11 & ageDeath < 20 ~ 12,
                              ageDeath > 20 ~ 24))

ref_geno <- "APOE4-KI_WT; Trem2-R47H_WT"

res_load1 <- get_all_de_results(
  meta_load1, counts,
  study = "Jax.IU.Pitt_APOE4.Trem2.R47H",
  ref_genotype = ref_geno,
  contrasts = list(
    c("genotype", "APOE4-KI_homozygous; Trem2-R47H_homozygous", ref_geno),
    c("genotype", "APOE4-KI_homozygous; Trem2-R47H_WT", ref_geno),
    c("genotype", "APOE4-KI_WT; Trem2-R47H_homozygous", ref_geno)),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("sequencingBatch", "real_age", "genotype")
)

# TODO write to file


## UCI_3xTg-AD -----------------------------------------------------------------

# There are no batch effects, no separate tissues
res_3x <- get_all_de_results(
  metadata_all, counts,
  study = "UCI_3xTg-AD",
  ref_genotype = "3xTg-AD_noncarrier",
  contrasts = list(c("genotype", "3xTg-AD_carrier", "3xTg-AD_noncarrier")),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("genotype")
)

# TODO write to file


## UCI_5XFAD -------------------------------------------------------------------

# This data has a sequencing batch, but it corresponds to the age of death so
# we don't include it in the model. There are two tissues.

res_5x <- get_all_de_results(
  metadata_all, counts, study = "UCI_5XFAD",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier")),
  group_cols = c("tissue", "sex", "ageDeath"),
  model_vars = c("genotype")
)
# TODO write to file


## UCI_ABCA7 -------------------------------------------------------------------

# TODO this data has rnaBatch with a bad split:
#   Abca7-V1599M_homozygous and 5XFAD_carrier; Abca7-V1599M_homozygous were
#     sequenced across two batches
#   5XFAD_carrier and 5XFAD_noncarrier were sequenced across two other batches
#   No overlap between the pairs of genotypes
# rnaBatch also has exact overlap with age of death, so each age of death contains
# one batch per pair of genotypes.
# Given this, we can't actually batch correct this data, which is not ideal.
res_abca7 <- get_all_de_results(
  metadata_all, counts,
  study = "UCI_ABCA7",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier"),
                   c("genotype", "Abca7-V1599M_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Abca7-V1599M_homozygous", "5XFAD_noncarrier")),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("genotype")
)
# TODO write to file


## UCI_hAbeta_KI ---------------------------------------------------------------

# There are no batches in this data
res_abki <- get_all_de_results(
  metadata_all, counts,
  study = "UCI_hAbeta_KI",
  ref_genotype = "hAbeta-KI_LoxP_WT",
  contrasts = list(c("genotype", "hAbeta-KI_LoxP_homozygous", "hAbeta-KI_LoxP_WT")),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("genotype")
)
# TODO write to file

## UCI_Trem2_Cuprizone ---------------------------------------------------------------

# This data only has 3-month males. No females or other ages.
# There are batches with this study (both libraryBatch and rnaBatch), though
# using rnaBatch is sufficient as all rnaBatches only contain one libraryBatch.
# All 5XFAD_noncarrier and Trem2-R47H_NSS_homozygous were run in the same batch,
# and all cuprizone Trem2-R47H_CSS_homozygous mice are also in that batch. A
# second batch contains some C57BL6J, some Trem2-KO, and some control
# Trem2-R47H_CSS_homozygous, while a third batch contains some C57BL6J, all
# control Trem2-KO, and some control Trem2-R47H_CSS_homozygous.
#
# It looks like the comparisons need to be:
#   Trem2-R47H_NSS_homozygous vs 5XFAD_noncarrier
#   cuprizone Trem2-R47H_CSS_homozygous vs cuprizone 5XFAD_noncarrier
#   Trem2-KO vs C57BL6J
#   control Trem2-R47H_CSS_homozygous vs control C57BL6J
#
# Or do we pool both C57BL6J and 5XFAD_noncarrier into one "WT" genotype?

res_trem2cup <- get_all_de_results(
  metadata_all, counts, study = "UCI_Trem2_Cuprizone",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list("??"),
  group_cols = c("treatmentType"),
  model_vars = c("rnaBatch", "genotype")
)


## UCI_Trem2-R47H_NSS ---------------------------------------------------------------

# TODO This data has all 3 batch variables filled in (libraryBatch,
# sequencingBatch, rnaBatch). None of them are entirely unique combinations of
# another even when split by age and sex. The best option is probably to paste
# all the batches together to create a new batch variable. When split by age,
# sex, and genotype, it looks like 5XFAD_noncarriers and 5XFAD_carriers were
# sequenced across 1-2 overlapping batches at each sex/time point, and
# Trem2-R47H_NSS_homozygous and 5XFAD_carrier; Trem2-R47H_NSS_homozygous were
# sequenced across 1-2 other overlapping batches, with no overlap between the
# two sets.
res_trem2nss <- get_all_de_results(
  metadata_all, counts,
  study = "UCI_Trem2-R47H_NSS",
  ref_genotype = "5XFAD_noncarrier",
  contrasts = list(c("genotype", "5XFAD_carrier", "5XFAD_noncarrier"),
                   c("genotype", "Trem2-R47H_NSS_homozygous", "5XFAD_noncarrier"),
                   c("genotype", "5XFAD_carrier; Trem2-R47H_NSS_homozygous", "5XFAD_noncarrier")),
  group_cols = c("sex", "ageDeath"),
  model_vars = c("batch", "genotype")
)



