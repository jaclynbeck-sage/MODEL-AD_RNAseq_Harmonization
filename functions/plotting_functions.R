library(dplyr)
library(ggplot2)
library(patchwork)

# Util functions ---------------------------------------------------------------

# Calculate a PCA of the data from a specific study using the 2000 most highly
# variable genes. Then, merge the first 10 PC values with the study's metadata
# and validation df to return a data frame that is easily plottable.
#
# Arguments:
#   validation_df - a data frame with one or both of "valid_variant" and
#     "valid_expression" booleans for each sample. This data frame should also
#     have metadata about the samples, including genotype.
#   counts - a matrix of counts data on the CPM scale
#   carrier_regex - a regex that will match the mutation "carrier" genotypes
#     when used with grepl
calculate_pca <- function(validation_df, counts, carrier_regex) {
  # Make sure both valid_variant and valid_expression are in the data frame
  if (!("valid_variant" %in% colnames(validation_df))) {
    validation_df$valid_variant <- TRUE
  }
  if (!("valid_expression" %in% colnames(validation_df))) {
    validation_df$valid_expression <- TRUE
  }

  # Convert counts to log2-scale
  counts_log <- log2(counts[, validation_df$unique_specimenID] + 0.5)

  # Use the 2000 most variable genes
  vars <- matrixStats::rowVars(as.matrix(counts_log))
  hv <- names(sort(vars, decreasing = TRUE))[1:2000]

  pc <- prcomp(t(counts_log[hv, ]))

  # Make a plottable data frame that incorporates the metadata
  pc_plot <- pc$x[, paste0("PC", 1:10)] |>
    merge(validation_df, by.x = "row.names", by.y = "unique_specimenID") |>
    dplyr::rename(unique_specimenID = Row.names) |>
    # Create some useful variables that can be used to color the plots
    mutate(
      has_mutation = grepl(carrier_regex, genotype),
      mismatch_type = case_when(
        valid_variant & valid_expression ~ "none",
        valid_variant & !valid_expression ~ "expression",
        !valid_variant & valid_expression ~ "variant",
        !valid_variant & !valid_expression ~ "expression + variant"
      )
    )

  return(pc_plot)
}


# Plotting ---------------------------------------------------------------------

# Given a data frame with PCA values and metadata, plot a grid of PCs vs
# covariates. Only (PC, covariate) pairs with > 0.5 correlation are plotted.
# NOTE: The plots may not be formatted nicely if there are long labels or too
# many plots.
# TODO consider doing 1 PC per column or per row for easier parsing
plot_pc_grid <- function(pc_plot, study_name, extra_vars = c()) {
  pc_sub <- pc_plot |>
    select(PC1:PC10, sex:sequencingBatch, -RIN,
           any_of(c("has_mutation", extra_vars))) |>
    mutate(ageDeath = paste(ageDeath, "months"))

  n_unique <- apply(pc_sub, 2, function(X) { length(na.omit(unique(X))) })
  pc_sub <- pc_sub[, n_unique >= 2]

  corr <- variancePartition::canCorPairs(
    paste("~", paste(colnames(pc_sub), collapse = " + ")),
    pc_sub,
    showWarnings = FALSE
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "cor") |>
    # Remove the diagonal
    subset(var1 != var2) |>
    # Keep anything with > 0.5 correlation
    subset(cor > 0.5)

  pairs <- lapply(which(grepl("PC", corr$var1)), function(row) {
    c(corr$var1[row], corr$var2[row])
  })

  plt_list <- lapply(pairs, function(pair) {
    pc <- pair[[1]]
    covariate <- pair[[2]]
    plt <- ggplot(pc_sub, aes(x = .data[[pc]], y = .data[[covariate]], color = .data[[covariate]])) +
      geom_jitter() +
      theme_bw() +
      theme(legend.position = "none")

    # force axis labels to be not aligned with each other across columns/rows
    free(plt, type = "label")
  })

  print(wrap_plots(plt_list) + plot_annotation(title = study_name))
}
