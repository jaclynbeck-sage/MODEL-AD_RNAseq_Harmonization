validate_3x <- function(metadata, geno_calls, counts, symbol_map,
                        genotype_pattern = "3xTg-AD_carrier") {
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls, genotype_pattern = genotype_pattern,
    gene_symbol_match = "APP", # MAPT not used in variant calling
    total_positions = 2) |>
    mutate(
      valid_3x_variant = (is_carrier & est_genotype == "carrier") |
        (!is_carrier & est_genotype != "carrier")
    )

  var_mismatches <- subset(valid_variants, !valid_3x_variant)
  print_variant_mismatches(var_mismatches, "3xTg-AD")

  # Check gene expression
  counts_df <- make_counts_df(metadata, counts, symbol_map,
                              c("APP", "MAPT")) %>%
    # Using > 1 CPM is sufficient for this genotype. For this study, all non-
    # carriers express 0 or extremely low levels of APP so we can use an "or"
    # operation here.
    mutate(expr_3x = APP > 1 | MAPT > 1)

  valid_expression <- subset(counts_df,
                             (grepl(genotype_pattern, genotype) & expr_3x) |
                               (!grepl(genotype_pattern, genotype) & !expr_3x))
  count_mismatches <- subset(counts_df,
                             !(merged_file_specimenID %in% valid_expression$merged_file_specimenID)) |>
    select(study_name, specimenID, merged_file_specimenID, genotype, expr_3x, APP, MAPT)

  print_expression_mismatches(count_mismatches, "3xTg-AD")

  counts_df$valid_3x_expression <- counts_df$merged_file_specimenID %in% valid_expression$merged_file_specimenID

  final <- merge(valid_variants, counts_df) |>
    mutate(valid = valid_3x_variant & valid_3x_expression)

  return(list(valid = select(final, merged_file_specimenID, valid),
              detail = select(final, -valid)))
}


validate_5x <- function(metadata, geno_calls, counts, symbol_map,
                        genotype_pattern = "5XFAD_carrier") {
  valid_variants <- get_variant_mismatches(
    metadata, geno_calls,
    genotype_pattern = genotype_pattern,
    gene_symbol_match = c("APP", "PSEN1"),
    total_positions = 6
  ) |>
    mutate(
      valid_5x_variant = (is_carrier & est_genotype == "carrier") |
             (!is_carrier & est_genotype != "carrier")
    )

  var_mismatches <- subset(valid_variants, !valid_5x_variant)
  print_variant_mismatches(var_mismatches, "5XFAD")

  # Also look at gene expression, because a lack of detected variant doesn't mean
  # the genotype isn't there.
  # For 5XFAD, expression of APP and PSEN1 in carriers should be really high, so
  # we consider variant calling to be pretty reliable in that all 6 detected
  # variants should be detected in carrier mice, while non-carrier mice should
  # have none but may have a few low-quality detections due to the homology
  # between the human and mouse genes. The following combinations of recorded and
  # estimated carrier status are valid:
  #   Recorded carrier genotype, estimated carrier genotype
  #   Recorded non-carrier genotype, estimated "ambiguous" or "non-carrier" genotype
  # Anything else (carrier + ambiguous/non-carrier, and non-carrier + carrier) is
  # not valid.
  counts_df <- make_counts_df(metadata, counts, symbol_map,
                              c("APP", "PSEN1")) |>
    # Using > 1 CPM is sufficient for this genotype. This statement needs to be an
    # "and" instead of an "or" due to the high level of homology between mouse App
    # and human APP leading to some (presumably) spurious counts getting mapped to
    # human APP in many WT samples, even those from studies unrelated to 5XFAD.
    mutate(expr_5x = APP > 1 & PSEN1 > 1)

  valid_expression <- subset(counts_df,
                             (grepl(genotype_pattern, genotype) & expr_5x) |
                               (!grepl(genotype_pattern, genotype) & !expr_5x))

  count_mismatches <- subset(counts_df,
                             !(merged_file_specimenID %in% valid_expression$merged_file_specimenID))

  print_expression_mismatches(count_mismatches, "5XFAD")

  counts_df$valid_5x_expression <- counts_df$merged_file_specimenID %in% valid_expression$merged_file_specimenID

  final <- merge(valid_variants, counts_df) |>
    mutate(valid = valid_5x_variant & valid_5x_expression)

  return(list(valid = select(final, merged_file_specimenID, valid),
              detail = select(final, -valid)))
}


validate_APOE4_KI <- function(metadata, counts, symbol_map,
                              genotype_pattern = "APOE4-KI_(homo|hetero)") {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "APOE") %>%
    # We need to use > 2 CPM as a threshold for this genotype: there are a few
    # non-carrier samples with > 1 but < 2 CPM expression that seem to match other
    # non-carriers and clearly don't match carriers in expression of mouse APOE,
    # so these are probably true non-carriers.
    mutate(expr_apoe = APOE > 2)

  valid_expression <- subset(counts_df,
                             (grepl(genotype_pattern, genotype) & expr_apoe) |
                               (!grepl(genotype_pattern, genotype) & !expr_apoe))

  count_mismatches <- subset(counts_df,
                             !(merged_file_specimenID %in% valid_expression$merged_file_specimenID))

  print_expression_mismatches(count_mismatches, "APOE4-KI")

  # TODO 3 of the 4 WT samples that fail from the LOAD2 study express APOE at ~20
  # CPM and should probably pass instead. On the other hand, these 3 look like
  # outliers on a PCA

  final <- counts_df |>
    mutate(valid_apoe4_expression = merged_file_specimenID %in% valid_expression$merged_file_specimenID,
           valid = valid_apoe4_expression)

  return(list(valid = select(final, merged_file_specimenID, valid),
              detail = select(final, -valid)))
}


validate_Abca7 <- function(metadata, geno_calls,
                           genotype_pattern = "Abca7-V1599M_(homo|hetero)") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Abca7-V1613M",
                                           total_positions = 3) |>
    mutate(valid = (is_carrier & est_genotype == "carrier") |
             (!is_carrier & est_genotype != "carrier"))

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "Abca7-V1599M")

  # There is no noticeable difference in Abca7 gene expression between genotypes
  # for each study so we do not validate based on expression.
  return(list(valid = select(valid_variants, merged_file_specimenID, valid),
              detail = valid_variants))
}


validate_Abi3 <- function(metadata, geno_calls,
                          genotype_pattern = "Abi3-S209F_homozygous") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Abi3-S212F",
                                           total_positions = 3) |>
    mutate(valid = is_carrier |
             (!is_carrier & est_genotype != "carrier"))

  # Only one carrier sample had a detected Abi3 variant. Given the relatively
  # low expression of Abi3 across samples, we do not assume that carrier
  # mismatches are actual mismatches since it's possible there just wasn't
  # enough coverage to hit the target region.

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "Abi3-S209F")

  # There is no noticeable difference in gene expression between genotypes so
  # we do not validate by expression.
  return(list(valid = select(valid_variants, merged_file_specimenID, valid),
              detail = valid_variants))
}


validate_Bin1 <- function() {
 # TODO
}


validate_CLU_KI <- function(metadata, counts, symbol_map,
                            genotype_pattern = "Clu-rs2279590_KI_homozygous") {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "CLU") |>
    # All but one non-carrier express 0 CLU, and the one non-carrier expresses
    # it at 0.12 CPM, but we use 1 just to be safe.
    mutate(expr_clu = CLU > 1)

  valid_expression <- subset(counts_df,
                             (grepl(genotype_pattern, genotype) & expr_clu) |
                               (!grepl(genotype_pattern, genotype) & !expr_clu))

  count_mismatches <- subset(counts_df,
                             !(merged_file_specimenID %in% valid_expression$merged_file_specimenID))

  print_expression_mismatches(count_mismatches, "CLU-KI")

  final <- counts_df |>
    mutate(valid = merged_file_specimenID %in% valid_expression$merged_file_specimenID)

  return(list(valid = select(final, merged_file_specimenID, valid),
              detail = final))
}


validate_hAbeta_KI <- function(metadata, geno_calls,
                               genotype_pattern = "hAbeta-KI_LoxP_homozygous") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "App-KI",
                                           total_positions = 3) |>
    mutate(valid = (is_carrier & est_genotype == "carrier")  |
             (!is_carrier & est_genotype != "carrier"))

  var_mismatches <- subset(valid_variants, !valid)
  print_variant_mismatches(var_mismatches, "hAbeta-KI")

  # There is no noticeable difference in App expression between genotypes, so we
  # do not validate by expression.
  return(list(valid = select(valid_variants, merged_file_specimenID, valid),
              detail = valid_variants))
}


validate_Picalm <- function() {
 # TODO
}


validate_Trem2_KO <- function(metadata, counts, symbol_map) {
  counts_df <- make_counts_df(metadata, counts, symbol_map, "Trem2") %>%
    # Threshold chosen by examination of plot of expression vs genotype
    mutate(expr_trem2 = Trem2 > 7)

  valid_expression <- subset(df_trem2ko, (genotype == "Trem2-KO" & !expr_trem2) |
                               (genotype != "Trem2-KO" & expr_trem2))
  count_mismatches <- subset(counts_df,
                             !(merged_file_specimenID %in% valid_expression$merged_file_specimenID))

  print_expression_mismatches(count_mismatches, "Trem2-KO")

  final <- counts_df |>
    mutate(valid = merged_file_specimenID %in% valid_expression$merged_file_specimenID)

  return(list(valid = select(final, merged_file_specimenID, valid),
              detail = final))
}


validate_Trem2_R47H <- function(metadata, geno_calls,
                                genotype_pattern = "Trem2-R47H(_NSS|_CSS)?_(homo|hetero)") {
  valid_variants <- get_variant_mismatches(metadata, geno_calls,
                                           genotype_pattern = genotype_pattern,
                                           mutation_match = "Trem2-R47H",
                                           total_positions = 1) |>
    mutate(valid_trem2_r47h_variant = is_carrier |
             (!is_carrier & est_genotype != "carrier"))

  # There is 1 carrier sample (homo- or heterozygous) with detected variants but
  # low quality, and another 5 carriers with no detected variants at all. The
  # breakdown is:
  #   low quality: 1 UCI_PrimaryScreen
  #   no detected variants: 2 Jax.IU.Pitt_APOE4.Trem2.R47H, 3 UCI_PrimaryScreen
  #
  # The 2 Jax samples have expression of Trem2 that is consistent with other
  # carrier samples, not with WT samples, so these 2 samples are probably valid.
  # Due to the lowered expression of Trem2 in R47H carriers, it's possible there
  # wasn't enough coverage to hit the target region enough, so we do not consider
  # carrier mismatches as true mismatches. Only non-carriers with detected
  # Trem2-R47H variants are considered mismatches.

  var_mismatches <- subset(valid_variants, !valid_trem2_r47h_variant)
  print_variant_mismatches(var_mismatches, "Trem2-R47H")

  # There is no clear difference in expression of Trem2 between genotypes in each
  # study so we can not validate based on expression.

  valid_variants <- mutate(valid_variants, valid = valid_trem2_r47h_variant)

  return(list(valid = select(valid_variants, merged_file_specimenID, valid),
              detail = select(valid_variants, -valid)))
}
