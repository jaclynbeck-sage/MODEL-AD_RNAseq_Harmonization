#### This code block contains libraries and general variables that should be
# loaded when this file is sourced
library(synapser)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(vsn)

# Synapse IDs used in the notebooks -- merged metadata and the gene_counts
# file with data from all studies
syn_metadata_file_id <- "syn61850266"
syn_gene_counts_id <- "syn62690577"
syn_symbol_map_id <- "syn62063692"

# Assumes this script is being sourced from inside the "visualization" folder
tmp_dir <- "tmp"

# ------------------------------------------------------------------------------


# Load the metadata and counts for a specific study, plus the mapping between
# Ensembl ID and gene symbol, and return them as a named list.
load_data <- function(study) {
  synLogin(silent = TRUE)

  metadata_file <- synGet(syn_metadata_file_id, downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")
  counts_file <- synGet(syn_gene_counts_id, downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")
  symbol_map_file <- synGet(syn_symbol_map_id, downloadLocation = tmp_dir,
                            ifcollision = "overwrite.local")

  metadata <- read.csv(metadata_file$path) %>%
    subset(study_name == study) %>%
    mutate(genotype = factor(genotype),
           sex = factor(sex),
           ageDeath = round(ageDeath, digits = 0))

  counts <- read.table(counts_file$path, header = TRUE, sep = "\t",
                       row.names = 1) %>%
    # Get rid of transcript_id column
    select(-transcript_id.s.) %>%
    # RSEM can return non-integer numbers for counts, so it needs to be rounded
    # to whole numbers for DESeq2
    round(digits = 0)

  # Subset to validated samples only
  # TODO temporary until this is finalized and on Synapse
  valid_samples <- read.csv("../data/Model_AD_valid_samples.csv")
  valid_samples <- subset(valid_samples, validated == TRUE)
  metadata <- subset(metadata,
                     merged_file_specimenID %in% valid_samples$merged_file_specimenID)

  # Not all samples in the metadata file appear in the counts matrix and vice versa
  metadata <- subset(metadata, merged_file_specimenID %in% colnames(counts))
  counts <- counts[, metadata$merged_file_specimenID]

  symbol_map <- read.csv(symbol_map_file$path)

  return(list("metadata" = metadata, "counts" = counts,
              "symbol_map" = symbol_map))
}


# Plot expression for the human transgenes
plot_human_genes <- function(counts_mat, metadata, symbol_map) {
  human_df <- counts_mat[grepl("ENSG", rownames(counts_mat)), ] %>%
    # Ensembl IDs disappear if we don't make it a column
    mutate(ensembl_gene_id = rownames(.)) %>%
    # All columns except ensembl_gene_id get melted into two columns, sample and
    # count
    tidyr::pivot_longer(cols = where(is.numeric), names_to = "sample",
                        values_to = "count") %>%
    # Get genotype information
    merge(metadata, by.x = "sample", by.y = "merged_file_specimenID") %>%
    # Get gene symbols
    merge(symbol_map)

  plt <- ggplot(human_df, aes(x = genotype, y = count, color = genotype)) +
    geom_boxplot(position = "dodge", outliers = FALSE) +
    geom_jitter(height = 0) +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~gene_symbol, scales = "free", nrow = 1) +
    scale_color_brewer(palette = "Set1")

  return(plt)
}
