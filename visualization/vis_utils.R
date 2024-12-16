#### This code block contains libraries and general variables that should be
# loaded when this file is sourced
library(synapser)
library(dplyr)
library(stringr)
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

# Load the metadata and counts for a specific study, convert Ensembl IDs to gene
# symbols, subset them to valid samples only, and return them as a named list.
# The list also includes a separate metadata and and a normalized counts df that
# includes samples that failed to pass validation.
load_data <- function(study) {
  synLogin(silent = TRUE)

  metadata_file <- synGet(syn_metadata_file_id,
                          downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

  counts_file <- synGet(syn_gene_counts_id,
                        downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")

  symbol_map_file <- synGet(syn_symbol_map_id,
                            downloadLocation = tmp_dir,
                            ifcollision = "overwrite.local")

  metadata <- read.csv(metadata_file$path) %>%
    subset(study_name == study) %>%
    mutate(genotype = factor(genotype),
           sex = factor(sex),
           ageDeath = round(ageDeath, digits = 0))

  counts <- read.table(counts_file$path, header = TRUE,
                       sep = "\t", row.names = 1) %>%
    # Get rid of transcript_id column
    select(-transcript_id.s.) %>%
    # RSEM can return non-integer numbers for counts, so it needs to be rounded
    # to whole numbers for DESeq2
    round(digits = 0)

  symbol_map <- read.csv(symbol_map_file$path)

  # Convert Ensembl IDs to gene symbols in the counts df
  rownames(symbol_map) <- symbol_map$ensembl_gene_id
  rownames(counts) <- make.unique(symbol_map[rownames(counts), "gene_symbol"])

  # Separate valid and invalid samples
  # TODO temporary until this file is finalized and on Synapse
  valid_samples <- read.csv("../data/Model_AD_valid_samples.csv")
  metadata <- merge(metadata, valid_samples) %>%
    dplyr::rename(valid_sample = validated)

  # Includes valid + invalid samples
  all_metadata <- metadata
  # Ensures sample order is the same for metadata and counts
  all_counts <- counts[, all_metadata$merged_file_specimenID]
  all_counts_norm <- vst(as.matrix(all_counts))

  # Valid samples only
  metadata <- subset(metadata, valid_sample == TRUE)
  counts <- counts[, metadata$merged_file_specimenID]

  return(list("metadata" = metadata,
              "counts" = counts,
              "all_metadata" = all_metadata,
              "all_counts" = all_counts,
              "all_counts_norm" = all_counts_norm))
}


# ------------------------------------------------------------------------------

print_invalid_samples <- function(all_metadata, columns_print) {
  if (sum(all_metadata$valid_sample) == nrow(all_metadata)) {
    print("All samples passed validation for this study.")
  } else {
    n_total_samples <- nrow(all_metadata)
    n_bad_samples <- n_total_samples - sum(all_metadata$valid_sample)
    print(
      str_glue(
        "{n_bad_samples} of {n_total_samples} samples failed validation ",
        "for this study."
      )
    )

    # Add standard columns used for all studies
    columns_print <- c("individualID", "specimenID", "sex", "genotype",
                       "valid_sex", columns_print)

    all_metadata %>%
      subset(!valid_sample) %>%
      select_at(columns_print)
  }
}

# ------------------------------------------------------------------------------

print_study_design <- function(metadata, group_split_vars, row_var, col_var) {
  if (length(group_split_vars) == 1) {
    metadata$group <- metadata[, group_split_vars]
  } else {
    metadata$group <- do.call(paste, metadata[, group_split_vars])
  }

  for (group_val in unique(metadata$group)) {
    meta_group <- subset(metadata, group == group_val)

    print(str_glue("{group_val}:"))
    print(table(meta_group[, row_var], meta_group[, col_var]))

    # Using str_glue makes a blank line without the [1] at the beginning
    # that usually prints out with print()
    print(str_glue(""))
  }
}

# ------------------------------------------------------------------------------

create_pca_plot <- function(metadata, counts_norm, title = "", n_hv = 500,
                            brewer_palette = "Set1", discrete_color = TRUE,
                            ...) {
  # X most highly-variable genes
  vars <- matrixStats::rowVars(counts_norm)
  hv <- names(sort(vars, decreasing = TRUE))[1:n_hv]

  pc <- prcomp(t(counts_norm[hv, ]))
  pc_plot <- merge(metadata, pc$x[, c("PC1", "PC2")],
                   by.x = "merged_file_specimenID",
                   by.y = "row.names")

  plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, ...)) +
    geom_point(size = 3) + ggtitle(title) + theme_bw()

  if (discrete_color) {
    plt <- plt + scale_color_brewer(palette = brewer_palette)
  } else {
    plt <- plt + scale_color_distiller(palette = brewer_palette, direction = 1)
  }

  plt
}

# ------------------------------------------------------------------------------

plot_gene_counts <- function(metadata, counts, genes_plot, draw_boxplot = TRUE,
                             draw_points = TRUE, title = "", ...) {
  plot_df <- counts[genes_plot, ] %>%
    as.data.frame() %>%
    mutate(gene = rownames(.)) %>%
    # All columns except "gene" get melted into two columns, sample and count
    tidyr::pivot_longer(cols = where(is.numeric),
                        names_to = "sample",
                        values_to = "counts") %>%
    # Get any metadata necessary for plotting
    merge(metadata, by.x = "sample", by.y = "merged_file_specimenID")

  plt <- ggplot(plot_df, aes(...))

  if (draw_boxplot) {
    plt <- plt + geom_boxplot(position = "dodge", outliers = FALSE)
  }
  if (draw_points) {
    plt <- plt + geom_jitter(height = 0)
  }

  plt <- plt + theme_bw() + ggtitle(title) +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap( ~ gene, scales = "free", nrow = 1) +
    scale_color_brewer(palette = "Set1")

  plt
}


plot_grouped_boxplot <- function(metadata, draw_boxplot = TRUE,
                                 draw_points = TRUE, title = "", ...) {
  plt <- ggplot(metadata, aes(...)) +
    ggtitle(title) + theme_bw() +
    scale_color_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (draw_boxplot) {
    plt <- plt + geom_boxplot(outliers = FALSE)
  }
  if (draw_points) {
    plt <- plt + geom_jitter(position = position_jitterdodge())
  }

  plt
}

# ------------------------------------------------------------------------------

shorten_genotypes <- function(geno_values, replacement_list) {
  for (pattern in names(replacement_list)) {
    geno_values <- str_replace_all(geno_values, pattern,
                                   replacement_list[[pattern]])
  }

  geno_values <- str_replace_all(geno_values, "; ", "/")
  return(geno_values)
}
