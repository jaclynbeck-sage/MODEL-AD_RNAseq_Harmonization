# Printing functions -----------------------------------------------------------

print_sex_mismatches <- function(xy_df) {
  mismatches <- subset(xy_df, !valid_sex) |>
    dplyr::rename(reported_sex = sex) |>
    mutate(across(c(Xist, Ddx3y, Eif2s3y, Kdm5d, mean_y), ~ round(.x, 2)))

  cat("*** Sex verification ***", "\n", "\n")

  cat(paste(nrow(mismatches), "samples have mismatched sex:"), "\n")
  print(select(mismatches, study, specimenID, reported_sex, est_sex,
               Xist, Eif2s3y, Ddx3y, Kdm5d, mean_y))

  cat("", "\n\n")
}

# Given a data frame of samples with genotype variant mismatches, print out this
# information grouped by carrier/non-carrier status.
print_variant_mismatches <- function(var_mismatches, genotype_name) {
  wt_mismatches <- subset(var_mismatches, !is_carrier) |> as.data.frame()
  carrier_mismatches <- subset(var_mismatches, is_carrier) |> as.data.frame()

  if (nrow(wt_mismatches) > 0) {
    cat(str_glue("Detected {genotype_name} variants in {nrow(wt_mismatches)} ",
                 "WT samples:"),
        "\n")

    print(select(wt_mismatches, study, specimenID, genotype, est_genotype))
  } else {
    cat(str_glue("No WT samples have detected {genotype_name} variants."),
        "\n")
  }

  cat("", "\n\n")

  if (nrow(carrier_mismatches) > 0) {
    cat(str_glue("{nrow(carrier_mismatches)} carrier samples had no detected ",
                 "{genotype_name} variants:"),
        "\n")
    print(select(carrier_mismatches, study, specimenID, genotype, est_genotype))
  } else {
    cat(str_glue("No carrier samples are missing detected {genotype_name} variants."),
        "\n")
  }

  cat("", "\n\n")
}


# Given a data frame of samples with expression mismatches, print out this
# information
# TODO group by carrier status?
print_expression_mismatches <- function(mismatch_df, genotype_name, genes) {
  if (nrow(mismatch_df) > 0) {
    cat(str_glue("{nrow(mismatch_df)} samples have a {genotype_name} expression mismatch:"),
        "\n")
    print(select(mismatch_df, -unique_specimenID))
  } else {
    cat(str_glue("No samples have a {genotype_name} expression mismatch."),
        "\n")
  }

  cat("", "\n\n")
}


# Utility functions ------------------------------------------------------------

# Bin quality scores from VCF files into "deletion", "low", "moderate", or
# "high" to help determine how reliable a detection is. Values of exactly 0
# were added by this code to indicate that no variant was detected at all.
score_quality <- function(quality) {
  case_when(
    is.na(quality) ~ "deletion",
    quality == 0 ~ "none",
    quality > 0 & quality < 20 ~ "low",
    quality >= 20 & quality < 100 ~ "moderate",
    quality >= 100 ~ "high"
  )
}


# Get a data frame with specific genes, merged with metadata. Adds one column
# per gene in "genes", with the column name being the gene symbol and the
# values being expression.
make_counts_df <- function(metadata, counts, symbol_map, genes) {
  genes_sub <- subset(symbol_map, gene_symbol %in% genes)
  rownames(genes_sub) <- genes_sub$ensembl_gene_id

  counts_sub <- counts[genes_sub$ensembl_gene_id, ]
  rownames(counts_sub) <- genes_sub[rownames(counts_sub), "gene_symbol"]

  counts_df <- as.data.frame(t(counts_sub)) |>
    merge(metadata, by.x = "row.names", by.y = "unique_specimenID") |>
    dplyr::rename(unique_specimenID = Row.names)
}


# Given a set of metadata and the variant detections for each sample, determine
# whether the genotype in the metadata matches the estimated genotype based on
# detected variants.
#
# Arguments:
#   metadata - data frame with one row per sample, which contains each mouse's
#     genotype. It is assumed that this data frame contains only the samples
#     we are validating and not all samples from all studies.
#   geno_info - a data frame with all variants detected in all samples. Samples
#     with no detected variants are not in this data frame
#   genotype_pattern - a regex used to grep for the "carrier" genotype, which
#     should only match genotypes that contain the modified gene and not any
#     wild-type genotypes.
#   mutation_match - optional, the name of the mutation(s) from the intervals.bed
#     file that should be included in this check. If NULL, the mutations are
#     selected using gene_symbol_match only. At least one of these two fields
#     should be non-NULL. mutation_match is typically used if there is only one
#     mutation for a given gene that is relevant to these samples (e.g. looking
#     at only the Trem2-R47H variant and not other variants for Trem2 in the
#     intervals file.)
#   gene_symbol_match - optional, the gene symbol of the mutation(s) that should
#     be included in this check. If NULL, the mutations are selected using
#     mutation_match only. At least one of these two fields should be non-NULL.
#     gene_symbol_match is typically used if there is more than one variant for
#     a gene and all variants should be included (e.g. all APP variants for
#     5xFAD mice).
#   total_positions - the number of expected genome positions where variants
#     should be detected in a "carrier" sample. This number should be the total
#     over all expected variants covered by mutation_match or gene_symbol_match.
#   rm_deletions - if the variants being judged are NOT deletions, any detected
#     deletions are probably false positives and are not counted toward the
#     number of variants detected in a sample.
get_variant_mismatches <- function(metadata, geno_info, genotype_pattern,
                                   mutation_match = NULL,
                                   gene_symbol_match = NULL,
                                   total_positions = 0,
                                   rm_deletions = TRUE) {

  # Subset to only variants for the specific gene or mutation specified, add
  # samples with no detected variants, then count the total number of detected
  # variants for each sample
  geno_sub <- geno_info |>
    subset(unique_specimenID %in% metadata$unique_specimenID) |>
    subset(gene_symbol %in% gene_symbol_match | mutation %in% mutation_match) |>

    # Change deletions to "none" if rm_deletions is TRUE, so they are not
    # counted below
    mutate(evidence = ifelse(rm_deletions & evidence == "deletion",
                             "none",
                             evidence)) |>

    # Summarize across all variants for each sample
    group_by(study, unique_specimenID, specimenID) |>

    # Count only variants with "moderate", "high", or (optionally) "deletion"
    # evidence / quality as detections. Keep track of average quality for debug
    # purposes.
    summarize(
      total_found = sum(evidence != "low" & evidence != "none"),
      avg_quality = mean(quality, na.rm = TRUE),
      evidence = paste(sort(unique(evidence)), collapse = ", "),
      .groups = "drop"
    ) |>

    # Merge with metadata to insert samples that don't exist in geno_sub
    # post-filtering (e.g. samples with no detected variants for this gene /
    # mutation), and to get the reported genotype
    merge(metadata, all.y = TRUE) |>

    # Determine genotype
    mutate(
      # Fill in NA evidence and total_found values from the merge
      evidence = ifelse(is.na(evidence), "none", evidence),
      total_found = ifelse(is.na(total_found), 0, total_found),

      # Mark which samples should have variants based on their reported genotype
      is_carrier = grepl(genotype_pattern, genotype),

      # Estimate genotype based on variants
      est_genotype = case_when(
        total_found >= total_positions ~ "carrier",
        total_found == 0 & evidence == "none" ~ "noncarrier",
        .default = "ambiguous"
      )
    )

  return(geno_sub)
}


# Merge two data frames that may have one or both of "valid_variant" and
# "valid_expression" columns. If they both have the same column, the merge creates
# new columns called "<col_name>.x" and "<col_name>.y". The values for .x and .y
# are combined with "&" so the final column is TRUE only if both .x and .y are TRUE.
# .x and .y columns are removed so there is only a single valid_variant and
# valid_expression column.
merge_validation_dfs <- function(df1, df2) {
  no_merge_columns <- c("valid_variant", "valid_expression", "is_carrier",
                        "est_genotype", "total_found", "avg_quality", "evidence")
  merged <- merge(df1, df2,
                  by = setdiff(intersect(colnames(df1), colnames(df2)),
                               no_merge_columns)) |>
    mutate(
      valid_variant = reduce(pick(starts_with("valid_variant")), `&`),
      valid_expression = reduce(pick(starts_with("valid_expression")), `&`)
    ) |>
    select(-ends_with(".x"), -ends_with(".y"))
}


# Sex validation functions -----------------------------------------------------

# Take the mean of several Y-chromosome genes and judge sex based on mean Y
# expression
validate_sex <- function(metadata, counts, symbol_map) {
  xy_df <- make_counts_df(metadata, counts, symbol_map,
                          c("Xist", "Eif2s3y", "Ddx3y", "Kdm5d")) |>
    # Use a small threshold for Xist for females -- there is one sample that
    # expresses near-zero counts of both Y-genes and Xist which should get
    # filtered out
    mutate(
      # Take the geometric mean (log-scale), including a pseudocount
      mean_y = rowMeans(log(cbind(Eif2s3y, Ddx3y, Kdm5d) + 0.5)),
      mean_y = exp(mean_y) - 0.5,
      # Expression of Y-related genes should be zero for females, so we use
      # anything > 1 CPM for males and assume < 1 CPM might be noise.
      est_sex = case_when(
        mean_y >= 1 ~ "male",
        Xist > 10 & mean_y < 1 ~ "female",
        .default = "unknown"
      ),
      # Mark valid/invalid sex matches
      valid_sex = (sex == est_sex)
    )

  print_sex_mismatches(xy_df)

  return(xy_df)
}


# Genotype validation functions ------------------------------------------------

## 3xTg-AD ---------------------------------------------------------------------

# Validate the presence or absence of APP-SweM671L / APP-SweK670N variants, and
# expression of human APP and MAPT in 3xTg-AD mice. The Psen1-M146V variant
# is not reliably detected in 3xTg-AD_carrier mice, so it is excluded from the
# genotype check.
validate_3x <- function(metadata, geno_calls, counts, symbol_map,
                        genotype_pattern = "3xTg-AD_carrier") {
  # APP variant detection
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls, genotype_pattern = genotype_pattern,
    gene_symbol_match = "APP",
    total_positions = 2
  ) |>
    mutate(
      valid_variant = (is_carrier & est_genotype == "carrier") |
        (!is_carrier & est_genotype != "carrier")
    )

  # Check gene expression of APP and MAPT
  counts_df <- make_counts_df(metadata, counts, symbol_map, c("APP", "MAPT")) |>
    # Using > 1 CPM is sufficient for this genotype. For this study, all non-
    # carriers express 0 or extremely low levels of APP and MAPT so we can use
    # an "or" operation here.
    mutate(expr_3x = APP > 1 | MAPT > 1,
           is_carrier = grepl(genotype_pattern, genotype),
           valid_expression = (is_carrier & expr_3x) | (!is_carrier & !expr_3x))

  # Combine the two results
  final <- merge(valid_variants, counts_df)

  # Print any mismatches
  print_variant_mismatches(subset(valid_variants, !valid_variant), "3xTg-AD")

  count_mismatches <- subset(counts_df, !valid_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_3x, APP, MAPT)
  print_expression_mismatches(count_mismatches, "3xTg-AD")

  return(final)
}


## 5XFAD -----------------------------------------------------------------------

# Validate the presence or absence of 6 FAD mutations (4 APP, 2 PSEN1) and
# expression of human APP and PSEN1 in 5XFAD mice.
#
# There are 6 total variant positions, but sometimes one of the PSEN1 positions
# has low-quality detection. We allow samples with at least 5 of the 6 variants
# to pass.
#
# Non-carrier mice should ideally have no detected variants, but some samples
# have a few low-quality detections due to the homology between the human and
# mouse genes. Non-carriers with an estimated genotype of "ambiguous" or
# "non-carrier" will pass.
validate_5x <- function(metadata, geno_calls, counts, symbol_map,
                        genotype_pattern = "5XFAD_carrier") {
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls,
    genotype_pattern = genotype_pattern,
    gene_symbol_match = c("APP", "PSEN1"),
    total_positions = 5 # pass with 5 of 6 positions
  ) |>
    mutate(
      valid_variant = (is_carrier & est_genotype == "carrier") |
             (!is_carrier & est_genotype != "carrier")
    )

  # Gene expression
  counts_df <- make_counts_df(metadata, counts, symbol_map, c("APP", "PSEN1")) |>
    # Using > 1 CPM is sufficient for this genotype. This statement needs to be an
    # "and" instead of an "or" due to the high level of homology between mouse App
    # and human APP leading to some (presumably) spurious counts getting mapped to
    # human APP in many WT samples, even those from studies unrelated to 5XFAD.
    mutate(expr_5x = APP > 1 & PSEN1 > 1,
           is_carrier = grepl(genotype_pattern, genotype),
           valid_expression = (is_carrier & expr_5x) | (!is_carrier & !expr_5x))

  # Combine the two data frames
  final <- merge(valid_variants, counts_df)

  # Print any mismatches
  print_variant_mismatches(subset(valid_variants, !valid_variant), "5XFAD")

  count_mismatches <- subset(counts_df, !valid_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_5x, APP, PSEN1)

  print_expression_mismatches(count_mismatches, "5XFAD")

  return(final)
}


## APOE4-KI --------------------------------------------------------------------

# Validate the presence or absence of human APOE expression in APOE4-KI mice.
# There are no detectable variants so this is judged by expression only.
#
# Unlike the other human genes, we need to use two different thresholds for this
# gene:
#   1. There is one LOAD2 sample that expresses 13 CPM of APOE and is a clear
#      outlier compared to APOE expression of other carriers, so we set a
#      threshold of > 20 CPM for positive expression of APOE.
#   2. There are a few non-carrier samples with > 1 but < 2 CPM expression
#      that seem to match other non-carriers and clearly don't match carriers
#      in expression of mouse Apoe, so these are probably true non-carriers.
#      However, there are also some non-carriers in the two LOAD2 studies that
#      express APOE at 7-20 CPM but express mouse Apoe at the same level as
#      other non-carriers. Three of these are also ambiguous sex mismatches,
#      so it's unclear whether these samples are contaminated or not. To be
#      cautious, we set the stricter limit of < 2 CPM for non-carriers.
validate_APOE4_KI <- function(metadata, counts, symbol_map,
                              genotype_pattern = "APOE4-KI_(homo|hetero)") {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "APOE") |>
    mutate(apoe_geno = grepl(genotype_pattern, genotype),
           valid_expression = (apoe_geno & APOE > 20) | (!apoe_geno & APOE < 2))

  # Print any mismatches
  count_mismatches <- subset(counts_df, !valid_expression) |>
    select(study, specimenID, unique_specimenID, genotype, APOE)

  print_expression_mismatches(count_mismatches, "APOE4-KI")

  return(counts_df)
}


## Abca7-V1599M ----------------------------------------------------------------

# Validate the presence or absence of the Abca7-V1599M variant in Abca7-V1599M
# mice. There is no noticeable difference in Abca7 gene expression between
# genotypes for each study so we do not validate based on expression.
validate_Abca7 <- function(metadata, geno_calls,
                           genotype_pattern = "Abca7-V1599M_(homo|hetero)") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Abca7-V1613M",
                                           total_positions = 3) |>
    mutate(valid_variant = (is_carrier & est_genotype == "carrier") |
             (!is_carrier & est_genotype != "carrier"))

  # Print any mismatches
  print_variant_mismatches(subset(valid_variants, !valid_variant), "Abca7-V1599M")

  return(valid_variants)
}


## Abi3-S209F ------------------------------------------------------------------

# Validate the presence or absence of the Abi3-S209F variant in Abi3-S209F mice.
# After testing, only one carrier sample had a detected Abi3 variant. Given the
# relatively low expression of Abi3 across samples, we do not assume that
# carrier mismatches are actual mismatches since it's possible there just wasn't
# enough coverage to hit the target region.
# There is no noticeable difference in gene expression between genotypes so
# we do not validate by expression.
validate_Abi3 <- function(metadata, geno_calls,
                          genotype_pattern = "Abi3-S209F_homozygous") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Abi3-S212F",
                                           total_positions = 3) |>
    mutate(valid_variant = is_carrier | # carriers pass even if variant not detected
             (!is_carrier & est_genotype != "carrier"))

  print_variant_mismatches(subset(valid_variants, !valid_variant), "Abi3-S209F")

  return(valid_variants)
}


## Bin1-K358R ------------------------------------------------------------------

# Validate the presence or absence of the Bin1-K358R variant in Bin1-K358R mice.
# There is no noticeable difference in expression between carriers and
# non-carriers, so this genotype is validated by variant detection only.
validate_Bin1 <- function(metadata, geno_calls, counts, symbol_map,
                          genotype_pattern = "Bin1-K358R_homozygous") {
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls,
    genotype_pattern = genotype_pattern,
    mutation_match = "Bin1-K358R",
    total_positions = 3
  ) |>
    mutate(
      valid_variant = (is_carrier & est_genotype == "carrier") |
        (!is_carrier & est_genotype != "carrier")
    )

  print_variant_mismatches(subset(valid_variants, !valid_variant), "UCI_Bin1K358R")

  return(valid_variants)
}


## Clu-rs2279590-KI ------------------------------------------------------------

# Validate the presence or absence of human CLU expression in Clu-rs2279590 mice.
# There are no detectable variants for this knock-in so we judge based on
# expression only.
validate_CLU_KI <- function(metadata, counts, symbol_map,
                            genotype_pattern = "Clu-rs2279590_KI_homozygous") {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "CLU") |>
    # All but one non-carrier express 0 CLU, and the one non-carrier expresses
    # it at 0.12 CPM, but we use 1 just to be safe.
    mutate(expr_clu = CLU > 1,
           is_carrier = grepl(genotype_pattern, genotype),
           valid_expression = (is_carrier & expr_clu) |
             (!is_carrier & !expr_clu))

  count_mismatches <- subset(counts_df, !valid_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_clu, CLU)

  print_expression_mismatches(count_mismatches, "CLU-KI")

  return(counts_df)
}


## hAPP-KI ---------------------------------------------------------------------

# Validate the presence or absence of the App knock-in variant in hAPP or
# hAbeta-KI mice. These mice do not express human APP, rather they have variants
# in mouse App. There is no noticeable difference in App expression between
# genotypes, so we do not validate by expression.
validate_hAPP_KI <- function(metadata, geno_calls,
                             genotype_pattern = "hAPP-3/3_homo|hetero") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "App-KI",
                                           total_positions = 3) |>
    mutate(valid_variant = (is_carrier & est_genotype == "carrier")  |
             (!is_carrier & est_genotype != "carrier"))

  print_variant_mismatches(subset(valid_variants, !valid_variant), "hAbeta-KI")

  return(valid_variants)
}


## Trem2-KO --------------------------------------------------------------------

# Validate the presence or absence of Trem2 KO by expression of Trem2.
validate_Trem2_KO <- function(metadata, counts, symbol_map) {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "Trem2") |>
    # Threshold chosen by examination of plot of expression vs genotype
    mutate(expr_trem2 = Trem2 > 7,
           valid_expression = (genotype == "Trem2-KO" & !expr_trem2) |
             (genotype != "Trem2-KO" & expr_trem2))

  count_mismatches <- subset(counts_df, !valid_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_trem2, Trem2)

  print_expression_mismatches(count_mismatches, "Trem2-KO")

  return(counts_df)
}


## Trem2-R47H ------------------------------------------------------------------

# Validate the presence or absence of Trem2-R47H variants in Trem2-R47H mice.
# This function works for both _NSS and _CSS variants. There is no clear cutoff
# in expression of Trem2 between genotypes in each study so we can not validate
# based on expression.
#
# Due to the generally lower expression of Trem2 in R47H carriers, it's possible
# there isn't enough coverage in carrier samples to hit the target region at
# detectable levels, so we do not consider lack of detected variants in carriers
# as true mismatches. Only non-carriers with detected Trem2-R47H variants are
# considered mismatches.
validate_Trem2_R47H <- function(metadata, geno_calls,
                                genotype_pattern = "Trem2-R47H(_NSS|_CSS)?_(homo|hetero)") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Trem2-R47H",
                                           total_positions = 1) |>
    mutate(valid_variant = is_carrier | # carriers always pass
             (!is_carrier & est_genotype != "carrier"))

  print_variant_mismatches(subset(valid_variants, !valid_variant), "Trem2-R47H")

  return(valid_variants)
}
