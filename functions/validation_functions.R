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
      valid_3x_variant = (is_carrier & est_genotype == "carrier") |
        (!is_carrier & est_genotype != "carrier")
    )

  # Check gene expression of APP and MAPT
  counts_df <- make_counts_df(metadata, counts, symbol_map, c("APP", "MAPT")) |>
    # Using > 1 CPM is sufficient for this genotype. For this study, all non-
    # carriers express 0 or extremely low levels of APP and MAPT so we can use
    # an "or" operation here.
    mutate(expr_3x = APP > 1 | MAPT > 1,
           is_carrier = grepl(genotype_pattern, genotype),
           valid_3x_expression = (is_carrier & expr_3x) | (!is_carrier & !expr_3x))

  # Combine the two results
  final <- merge(valid_variants, counts_df) |>
    mutate(valid = valid_3x_variant & valid_3x_expression)

  # Print any mismatches
  var_mismatches <- subset(valid_variants, !valid_3x_variant)
  print_variant_mismatches(var_mismatches, "3xTg-AD")

  count_mismatches <- subset(counts_df, !valid_3x_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_3x, APP, MAPT)
  print_expression_mismatches(count_mismatches, "3xTg-AD")

  return(list(valid = select(final, unique_specimenID, valid),
              detail = select(final, -valid)))
}


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
      valid_5x_variant = (is_carrier & est_genotype == "carrier") |
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
           valid_5x_expression = (is_carrier & expr_5x) | (!is_carrier & !expr_5x))

  # Combine the two data frames
  final <- merge(valid_variants, counts_df) |>
    mutate(valid = valid_5x_variant & valid_5x_expression)

  # Print any mismatches
  var_mismatches <- subset(valid_variants, !valid_5x_variant)
  print_variant_mismatches(var_mismatches, "5XFAD")

  count_mismatches <- subset(counts_df, !valid_5x_expression) |>
    select(study, specimenID, unique_specimenID, genotype, expr_5x, APP, PSEN1)

  print_expression_mismatches(count_mismatches, "5XFAD")

  return(list(valid = select(final, unique_specimenID, valid),
              detail = select(final, -valid)))
}


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
           valid_apoe4_expression = (apoe_geno & APOE > 20) | (!apoe_geno & APOE < 2),
           valid = valid_apoe4_expression)

  # Print any mismatches
  count_mismatches <- subset(counts_df, !valid_apoe4_expression) |>
    select(study, specimenID, unique_specimenID, genotype, APOE)

  print_expression_mismatches(count_mismatches, "APOE4-KI")

  return(list(valid = select(counts_df, unique_specimenID, valid),
              detail = select(counts_df, -valid)))
}


# Validate the presence or absence of the Abca7-V1599M variant in Abca7-V1599M
# mice. There is no noticeable difference in Abca7 gene expression between
# genotypes for each study so we do not validate based on expression.
validate_Abca7 <- function(metadata, geno_calls,
                           genotype_pattern = "Abca7-V1599M_(homo|hetero)") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Abca7-V1613M",
                                           total_positions = 3) |>
    mutate(valid = (is_carrier & est_genotype == "carrier") |
             (!is_carrier & est_genotype != "carrier"))

  # Print any mismatches
  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "Abca7-V1599M")

  return(list(valid = select(valid_variants, unique_specimenID, valid),
              detail = valid_variants))
}


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
    mutate(valid = is_carrier | # carriers pass even if variant not detected
             (!is_carrier & est_genotype != "carrier"))

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "Abi3-S209F")

  return(list(valid = select(valid_variants, unique_specimenID, valid),
              detail = valid_variants))
}


# Validate the presence or absence of the Bin1-K358R variant in Bin1-K358R mice.
validate_Bin1 <- function(metadata, geno_calls, counts, symbol_map,
                          genotype_pattern = "Bin1-K358R_homozygous") {
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls,
    genotype_pattern = genotype_pattern,
    mutation_match = "Bin1-K358R",
    total_positions = 3
  ) |>
    mutate(
      valid = (is_carrier & est_genotype == "carrier") |
        (!is_carrier & est_genotype != "carrier")
    )

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "UCI_Bin1K358R")

  return(list(valid = select(valid_variants, unique_specimenID, valid),
              detail = valid_variants))
}


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
           valid = (is_carrier & expr_clu) |
             (!is_carrier & !expr_clu))

  count_mismatches <- subset(counts_df, !valid) |>
    select(study, specimenID, unique_specimenID, genotype, expr_clu, CLU)

  print_expression_mismatches(count_mismatches, "CLU-KI")

  return(list(valid = select(counts_df, unique_specimenID, valid),
              detail = counts_df))
}


# Validate the presence or absence of the App knock-in variant in hAPP or
# hAbeta-KI mice. These mice do not express human APP, rather they have variants
# in mouse App. There is no noticeable difference in App expression between
# genotypes, so we do not validate by expression.
validate_hAbeta_KI <- function(metadata, geno_calls,
                               genotype_pattern = "hAbeta-KI_LoxP_homozygous") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "App-KI",
                                           total_positions = 3) |>
    mutate(valid_hAbeta_variant = (is_carrier & est_genotype == "carrier")  |
             (!is_carrier & est_genotype != "carrier"),
           valid = valid_hAbeta_variant)

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "hAbeta-KI")

  return(list(valid = select(valid_variants, unique_specimenID, valid),
              detail = select(valid_variants, -valid)))
}


validate_Picalm <- function() {
 # TODO
}


# Validate the presence or absence of Trem2 KO by expression of Trem2.
validate_Trem2_KO <- function(metadata, counts, symbol_map) {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "Trem2") |>
    # Threshold chosen by examination of plot of expression vs genotype
    mutate(expr_trem2 = Trem2 > 7,
           valid = (genotype == "Trem2-KO" & !expr_trem2) |
             (genotype != "Trem2-KO" & expr_trem2))

  count_mismatches <- subset(counts_df, !valid) |>
    select(study, specimenID, unique_specimenID, genotype, expr_trem2, Trem2)

  print_expression_mismatches(count_mismatches, "Trem2-KO")

  return(list(valid = select(counts_df, unique_specimenID, valid),
              detail = counts_df))
}


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
    mutate(valid_trem2_r47h_variant = is_carrier | # carriers always pass
             (!is_carrier & est_genotype != "carrier"),
           valid = valid_trem2_r47h_variant)

  var_mismatches <- subset(valid_variants, !valid_trem2_r47h_variant)
  print_variant_mismatches(var_mismatches, "Trem2-R47H")

  return(list(valid = select(valid_variants, unique_specimenID, valid),
              detail = select(valid_variants, -valid)))
}
