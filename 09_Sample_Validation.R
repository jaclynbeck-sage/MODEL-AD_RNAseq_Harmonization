library(synapser)
library(synapserutils)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

# Set up -----------------------------------------------------------------------

# Synapse IDs used in this script -- merged metadata and the gene counts file
# with data from all studies
syn_metadata_file_id <- "syn61850266"
syn_gene_counts_id <- "syn62690577"
syn_symbol_map_id <- "syn62063692"
syn_genotype_folder_id <- "syn63913842"
syn_intervals_id <- "syn69046595"

synLogin(silent = TRUE)

github <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/06_Sample_Validation.Rmd"
tmp_dir <- file.path("data", "tmp")


# Utility functions ------------------------------------------------------------

# Get a data frame with specific genes, merged with metadata
make_counts_df <- function(metadata, counts, symbol_map, genes) {
  genes_sub <- subset(symbol_map, gene_symbol %in% genes)
  rownames(genes_sub) <- genes_sub$ensembl_gene_id

  counts_sub <- counts[genes_sub$ensembl_gene_id, ]
  rownames(counts_sub) <- genes_sub[rownames(counts_sub), "gene_symbol"]

  counts_df <- as.data.frame(t(counts_sub)) %>%
    merge(metadata, by.x = "row.names", by.y = "merged_file_specimenID") %>%
    dplyr::rename(merged_file_specimenID = Row.names)
}


score_quality <- function(quality) {
  ifelse(is.na(quality), "deletion",
         ifelse(quality < 20, "low",
                ifelse(quality < 100, "moderate",
                       "high")))
}


get_variant_mismatches <- function(metadata, geno_info, genotype_pattern,
                                   mutation_match = NULL,
                                   gene_symbol_match = NULL,
                                   total_positions = 0,
                                   rm_na = TRUE) {
  carrier_samples <- subset(metadata, grepl(genotype_pattern, genotype))

  # Get every sample in the same studies as the carrier samples
  study_samples <- subset(metadata, study_name %in% carrier_samples$study_name)

  geno_sub <- subset(geno_info, merged_file_specimenID %in% study_samples$merged_file_specimenID &
                       (gene_symbol %in% gene_symbol_match |
                          mutation %in% mutation_match)) %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype, gene_symbol) %>%
    summarize(count = n(),
              avg_quality = mean(quality, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(evidence = score_quality(avg_quality),
           is_carrier = grepl(genotype_pattern, genotype))

  # Add in samples that don't show up in geno_sub
  missing_samples <- study_samples %>%
    subset(!(merged_file_specimenID %in% geno_sub$merged_file_specimenID)) %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype) %>%
    reframe(gene_symbol = unique(geno_sub$gene_symbol),
            count = 0, avg_quality = 0, evidence = "none",
            is_carrier = grepl(genotype_pattern, genotype))

  geno_sub <- rbind(geno_sub, missing_samples)

  geno_final <- geno_sub %>%
    group_by(study_name, merged_file_specimenID, specimenID, genotype, is_carrier) %>%
    summarize(total_found = sum(count),
              evidence = paste(unique(evidence), collapse = ","),
              .groups = "drop")

  if (rm_na) {
    geno_final$total_found[geno_final$evidence == "deletion"] <- 0
  }

  geno_final <- geno_final %>%
    mutate(est_genotype =
             ifelse(total_found == total_positions, "carrier",
                    ifelse(total_found == 0 & evidence == "none", "noncarrier",
                           "ambiguous")))

  return(geno_final)
}


print_variant_mismatches <- function(var_mismatches, genotype_name) {
  wt_mismatches <- subset(var_mismatches, !is_carrier)
  carrier_mismatches <- subset(var_mismatches, is_carrier)

  if (nrow(wt_mismatches) > 0) {
    print(paste("Detected", genotype_name, "variants in",
                nrow(wt_mismatches), "WT samples:"))

    print(select(wt_mismatches, study_name, specimenID, genotype, est_genotype))
  } else {
    print(paste("No WT samples have detected", genotype_name, "variants."))
  }

  print("")

  if (nrow(carrier_mismatches) > 0) {
    print(paste(nrow(carrier_mismatches), "carrier samples",
                "had no detected", genotype_name, "variants:"))
    print(select(carrier_mismatches, study_name, specimenID, genotype, est_genotype))
  } else {
    print(paste("No carrier samples are missing detected", genotype_name,
                "variants."))
  }
}


print_expression_mismatches <- function(mismatch_df, genotype_name) {
  if (nrow(mismatch_df) > 0) {
    print(paste(nrow(mismatch_df), "samples have a", genotype_name,
                "expression mismatch:"))
    print(select(mismatch_df, -merged_file_specimenID))
  } else {
    print(paste("No samples have a", genotype_name, "expression mismatch."))
  }
}


# Load counts and metadata -----------------------------------------------------

metadata_file <- synGet(syn_metadata_file_id,
                        downloadLocation = tmp_dir,
                        ifcollision = "overwrite.local")
counts_file <- synGet(syn_gene_counts_id,
                      downloadLocation = tmp_dir,
                      ifcollision = "overwrite.local")
symbol_map_file <- synGet(syn_symbol_map_id,
                          downloadLocation = tmp_dir,
                          ifcollision = "overwrite.local")

metadata_all <- read.csv(metadata_file$path)
counts <- read.delim(counts_file$path, header = TRUE, row.names = 1) %>%
  # Get rid of transcript_id column
  select(-transcript_id.s.)

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, merged_file_specimenID %in% colnames(counts)) %>%
  select(individualID, specimenID, merged_file_specimenID, sex, genotype, ageDeath, study_name)
counts <- counts[, metadata_all$merged_file_specimenID]

symbol_map <- read.csv(symbol_map_file$path)

# Convert counts to CPM so samples are comparable to each other
lib_sizes <- colSums(counts)
counts <- sweep(counts, 2, lib_sizes, "/") * 1e6


# Find vcf files on Synapse ----------------------------------------------------

# TODO: ideally these would be annotated with the specimenID / study name but
# for now we walk the folder structure and use the folder name for each sample
# to determine the specimenID

studies <- synGetChildren(syn_genotype_folder_id, includeTypes = list("folder"))$asList()
names(studies) <- sapply(studies, "[[", "name")

studies <- studies[names(studies) %in% metadata_all$study_name]

study_vcfs <- lapply(studies, function(study) {
  folder_structure <- synapserutils::walk(study$id, includeTypes = list("file"))$asList()

  samples_df <- lapply(folder_structure, function(item) {
    # The structure of each item that contains an actual file is a list:
    #   item[[1]] = list(1 = path of containing folder, 2 = synapse id of folder)
    #   item[[2]] = list(empty), or a list of sub-folders
    #   item[[3]] = list of files, where each item is list(1 = file name, 2 = synapse id)
    if (length(item) < 3 || lengths(item)[3] == 0) {
      return(NULL)
    }

    folder_name <- pluck(item, 1, 1)
    specimenID <- basename(folder_name)

    # Find the file that ends in .gz, not .gz.tbi
    vcf_filenames <- sapply(item[[3]], "[[", 1)
    ind <- which(grepl("gz$", vcf_filenames))
    return(data.frame(study_name = study$name,
                      specimenID = specimenID,
                      filename = vcf_filenames[ind],
                      syn_id = pluck(item, 3, ind, 2)))
  })

  samples_df <- do.call(rbind, samples_df)
  return(samples_df)
})


# Download and parse vcf files -------------------------------------------------

# Read the intervals file used to call variants and turn it into one row per
# genome position rather than one row per range of positions.
# The intervals file coordinates are -1 and +1 off from the real coordinate(s)
# as reported in the vcf file:
#   start = real coordinate - 1 and
#   end = real coordinate + 1
# We undo this operation so the values in the "position" column match the values
# in the vcf files.
intervals_file <- synGet(syn_intervals_id,
                         downloadLocation = tmp_dir,
                         ifcollision = "overwrite.local")

intervals <- read.delim(intervals_file$path, header = FALSE) %>%
  dplyr::rename(chromosome = V1, start = V2, end = V3, mutation = V4) %>%
  rowwise() %>%
  # adds as many new rows as necessary to cover the full range of positions
  reframe(chromosome = chromosome,
          mutation = mutation,
          position = (start+1):(end-1))

geno_info <- lapply(study_vcfs, function(study_df) {
  mutation_df <- apply(study_df, 1, function(row) {
    vcf_file <- synGet(row[["syn_id"]], downloadLocation = tmp_dir,
                       ifcollision = "overwrite.local")

    vcf <- NULL

    # If the sample didn't have any detected mutations, the vcf file will
    # contain only comments and no data, which will cause the read.delim
    # function to throw an error.
    try({
      vcf <- read.delim(gzfile(vcf_file$path), header = FALSE, comment.char = "#") %>%
        dplyr::rename(chromosome = V1, position = V2, ref = V4, alt = V5,
                      quality = V6) %>%
        select(chromosome, position, ref, alt, quality) %>%
        mutate(study_name = row[["study_name"]],
               specimenID = row[["specimenID"]]) %>%
        merge(intervals)
    }, silent = TRUE)

    # If no mutations found, this will be NULL
    return(vcf)
  })

  mutation_df <- do.call(rbind, mutation_df) %>%
    mutate(gene_symbol = str_replace(mutation, "-.*", ""))

  return(mutation_df)
})

geno_info <- do.call(rbind, geno_info)
geno_info$quality[geno_info$quality == "."] <- NA
geno_info$quality <- as.numeric(geno_info$quality)
geno_info <- merge(geno_info, metadata_all) %>%
  mutate(evidence = score_quality(quality))


# Set up the valid samples data frame ------------------------------------------

valid_samples <- metadata_all %>%
  select(individualID, specimenID, merged_file_specimenID, study_name)


# Validation of mouse sex ------------------------------------------------------

xy_df <- make_counts_df(metadata_all, counts, symbol_map,
                        c("Xist", "Eif2s3y", "Ddx3y")) %>%
  # More than 10 CPM for Xist -- there is low-level expression even in males
  mutate(est_female = Xist > 10,
         # Expression of Y-related genes should be zero for females, so we
         # use anything > 1 CPM for males and assume < 1 CPM might be noise
         est_male = Eif2s3y > 1 | Ddx3y > 1)

# Mark valid/invalid conditions
xy_df$valid_sex <- FALSE
xy_df$valid_sex[xy_df$sex == "female" &
                  xy_df$est_female == TRUE &
                  xy_df$est_male == FALSE] <- TRUE
xy_df$valid_sex[xy_df$sex == "male" &
                  xy_df$est_female == FALSE &
                  xy_df$est_male == TRUE] <- TRUE

# Mark valid samples in the data frame
sex_matches <- select(xy_df, merged_file_specimenID, valid_sex)
valid_samples <- merge(valid_samples, sex_matches, all = TRUE)

# For printing
mismatches <- subset(xy_df, !valid_sex) %>%
  select(study_name, individualID, specimenID, sex, est_female, est_male, Xist, Eif2s3y, Ddx3y) %>%
  dplyr::rename(reported_sex = sex)
mismatches$est_sex <- "unknown"
mismatches$est_sex[mismatches$est_female & !mismatches$est_male] <- "female"
mismatches$est_sex[mismatches$est_male & !mismatches$est_female] <- "male"

print(paste(nrow(mismatches), "samples have mismatched sex:"))
print(select(mismatches, study_name, specimenID, reported_sex, est_sex,
             Xist, Eif2s3y, Ddx3y))


# Validation of mouse genotype -------------------------------------------------

## 5XFAD carrier / non-carriers ------------------------------------------------

variants_5x <- get_variant_mismatches(metadata_all, geno_info,
                                      genotype_pattern = "5XFAD_carrier",
                                      gene_symbol_match = c("APP", "PSEN1"),
                                      total_positions = 6)

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
variants_merge_5x <- variants_5x %>%
  mutate(valid_5x_variant = (is_carrier & est_genotype == "carrier") |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_5x_variant)

valid_samples <- merge(valid_samples, variants_merge_5x, all = TRUE)

# Samples that weren't from 5X-related studies will have NA values in this
# column, so they should be valid by default.
valid_samples$valid_5x_variant[is.na(valid_samples$valid_5x_variant)] <- TRUE

var_mismatches_5x <- merge(variants_5x, subset(valid_samples, !valid_5x_variant))
print_variant_mismatches(var_mismatches_5x, "5XFAD")

# Also look at gene expression, because a lack of detected variant doesn't mean
# the genotype isn't there
studies_5x <- subset(metadata_all, grepl("5XFAD_carrier", genotype))
meta_5x <- subset(metadata_all, study_name %in% studies_5x$study_name)

df_5x <- make_counts_df(meta_5x, counts, symbol_map,
                        c("APP", "PSEN1")) %>%
  # Using > 1 CPM is sufficient for this genotype. This statement needs to be an
  # "and" instead of an "or" due to the high level of homology between mouse App
  # and human APP leading to some (presumably) spurious counts getting mapped to
  # human APP in many WT samples, even those from studies unrelated to 5XFAD.
  mutate(expr_5x = APP > 1 & PSEN1 > 1)

valid_5x <- subset(df_5x, (grepl("5XFAD_carrier", genotype) & expr_5x) |
                     (!grepl("5XFAD_carrier", genotype) & !expr_5x))
mismatch_5x <- subset(df_5x, !(merged_file_specimenID %in% valid_5x$merged_file_specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, expr_5x, APP, PSEN1)

print_expression_mismatches(mismatch_5x, "5XFAD")

valid_samples$valid_5x_expression <- !(valid_samples$merged_file_specimenID %in%
                                         mismatch_5x$merged_file_specimenID)

# Tmp PCA

counts_5x <- log2(counts[, meta_5x$merged_file_specimenID] + 0.5)

for (study in unique(meta_5x$study_name)) {
  meta_study <- subset(meta_5x, study_name == study)
  vars <- matrixStats::rowVars(as.matrix(counts_5x[, meta_study$merged_file_specimenID]))
  hv <- names(sort(vars, decreasing = TRUE))[1:2000]
  pc <- prcomp(t(counts_5x[hv, meta_study$merged_file_specimenID]))
  pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                   by.y = "merged_file_specimenID") %>%
    dplyr::rename(merged_file_specimenID = Row.names) %>%
    mutate(has_5x = grepl("5XFAD_carrier", genotype),
           mismatch_type = ifelse(merged_file_specimenID %in% mismatch_5x$merged_file_specimenID, "expression",
                                  ifelse(merged_file_specimenID %in% var_mismatches_5x$merged_file_specimenID, "variant",
                                         "none")))
  plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = has_5x)) +
    geom_point() + ggtitle(study)
  print(plt)
}


## 3xTg carrier / non-carriers -------------------------------------------------

variants_3x <- get_variant_mismatches(metadata_all, geno_info,
                                      genotype_pattern = "3xTg-AD_carrier",
                                      gene_symbol_match = "APP", # MAPT not used in variant calling
                                      total_positions = 2)

# Same concept as 5XFAD validation
variants_merge_3x <- variants_3x %>%
  mutate(valid_3x_variant = (is_carrier & est_genotype == "carrier") |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_3x_variant)

valid_samples <- merge(valid_samples, variants_merge_3x, all = TRUE)
valid_samples$valid_3x_variant[is.na(valid_samples$valid_3x_variant)] <- TRUE

var_mismatches_3x <- merge(variants_3x, subset(valid_samples, !valid_3x_variant))
print_variant_mismatches(var_mismatches_3x, "3xTg-AD")

# Check gene expression
studies_3x <- subset(metadata_all, grepl("3xTg-AD_carrier", genotype))
meta_3x <- subset(metadata_all, study_name %in% studies_3x$study_name)
df_3x <- make_counts_df(meta_3x, counts, symbol_map,
                        c("APP", "MAPT")) %>%
  # Using > 1 CPM is sufficient for this genotype. For this study, all non-
  # carriers express 0 or extremely low levels of APP so we can use an "or"
  # operation here.
  mutate(expr_3x = APP > 1 | MAPT > 1)

valid_3x <- subset(df_3x, (grepl("3xTg-AD_carrier", genotype) & expr_3x) |
                     (!grepl("3xTg-AD_carrier", genotype) & !expr_3x))
mismatch_3x <- subset(df_3x, !(merged_file_specimenID %in% valid_3x$merged_file_specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, expr_3x, APP, MAPT)

print_expression_mismatches(mismatch_3x, "3xTg-AD")

valid_samples$valid_3x_expression <- !(valid_samples$merged_file_specimenID %in%
                                         mismatch_3x$merged_file_specimenID)


## APOE4-KI homozygous / noncarrier --------------------------------------------

# We can't detect APOE4 KI with variant calling, just with gene expression. The
# Trem2 variants will be detected separately
studies_apoe <- subset(metadata_all, grepl("APOE4-KI", genotype))
meta_apoe <- subset(metadata_all, study_name %in% studies_apoe$study_name)
df_apoe <- make_counts_df(meta_apoe, counts, symbol_map, "APOE") %>%
  # We need to use > 2 CPM as a threshold for this genotype: there are a few
  # non-carrier samples with > 1 but < 2 CPM expression that seem to match other
  # non-carriers and clearly don't match carriers in expression of mouse APOE,
  # so these are probably true non-carriers.
  mutate(expr_apoe = APOE > 2)

valid_apoe <- subset(df_apoe, (grepl("APOE4-KI_(homo|hetero)", genotype) & expr_apoe) |
                     (!grepl("APOE4-KI_(homo|hetero)", genotype) & !expr_apoe))
mismatch_apoe <- subset(df_apoe, !(merged_file_specimenID %in% valid_apoe$merged_file_specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, expr_apoe, APOE)

print_expression_mismatches(mismatch_apoe, "APOE4-KI")

valid_samples$valid_apoe <- !(valid_samples$merged_file_specimenID %in%
                                mismatch_apoe$merged_file_specimenID)

# Tmp PCA

counts_apoe <- log2(counts[, meta_apoe$merged_file_specimenID] + 0.5)

for (study in unique(meta_apoe$study_name)) {
  meta_study <- subset(meta_apoe, study_name == study & !is.na(genotype))
  vars <- matrixStats::rowVars(as.matrix(counts_apoe[, meta_study$merged_file_specimenID]))
  hv <- names(sort(vars, decreasing = TRUE))[1:2000]
  pc <- prcomp(t(counts_apoe[hv, meta_study$merged_file_specimenID]))
  pc_plot <- merge(pc$x[, c("PC1", "PC2")], meta_study, by.x = "row.names",
                   by.y = "merged_file_specimenID") %>%
    dplyr::rename(merged_file_specimenID = Row.names) %>%
    mutate(has_apoe = grepl("APOE4-KI_(homo|hetero)", genotype),
           mismatch_type = ifelse(merged_file_specimenID %in% mismatch_apoe$merged_file_specimenID, "expression",
                                  "none"))
  plt <- ggplot(pc_plot, aes(x = PC1, y = PC2, color = mismatch_type, shape = has_apoe)) +
    geom_point() + ggtitle(study)
  print(plt)
}


## Abca7-V1599M ----------------------------------------------------------------

variants_abca7 <- get_variant_mismatches(metadata_all, geno_info,
                                         genotype_pattern = "Abca7-V1599M_(homo|hetero)",
                                         mutation_match = "Abca7-V1613M",
                                         total_positions = 3)

# We currently don't have any samples with an estimated genotype of "ambiguous".
# The 2 carrier mismatches correspond to 2 of the 3 carrier samples in the
# UCI_PrimaryScreen study, with no detected variants. The UCI_PrimaryScreen
# study seems to have less sequencing depth than other studies, and the
# expression of Abca7 is much lower across all samples in that study than in the
# UCI_ABCA7 study, so it's possible that there wasn't enough coverage to hit the
# target region in this sample. The UCI_PrimaryScreen study is excluded from
# further analysis so we will keep stringent validation criteria and re-evaluate
# if another study with Abca7 mutants gets added.
variants_merge_abca7 <- variants_abca7 %>%
  mutate(valid_abca7_variant = (is_carrier & est_genotype == "carrier") |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_abca7_variant)

valid_samples <- merge(valid_samples, variants_merge_abca7, all = TRUE)
valid_samples$valid_abca7_variant[is.na(valid_samples$valid_abca7_variant)] <- TRUE

var_mismatches_abca7 <- merge(variants_abca7, subset(valid_samples, !valid_abca7_variant))
print_variant_mismatches(var_mismatches_abca7, "Abca7-V1599M")

# There is no noticeable difference in Abca7 gene expression between genotypes
# for each study so we do not validate based on expression.


## Abi3-S209F ------------------------------------------------------------------

variants_abi3 <- get_variant_mismatches(metadata_all, geno_info,
                                        genotype_pattern = "Abi3-S209F_homozygous",
                                        mutation_match = "Abi3-S212F",
                                        total_positions = 3)

# Only one carrier sample had a detected Abi3 variant. Given the relatively low
# expression of Abi3 across samples, so we do not assume that carrier mismatches
# are actual mismatches since it's possible there just wasn't enough coverage to
# hit the target region.
variants_merge_abi3 <- variants_abi3 %>%
  mutate(valid_abi3_variant = is_carrier  |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_abi3_variant)

valid_samples <- merge(valid_samples, variants_merge_abi3, all = TRUE)
valid_samples$valid_abi3_variant[is.na(valid_samples$valid_abi3_variant)] <- TRUE

var_mismatches_abi3 <- merge(variants_abi3, subset(valid_samples, !valid_abi3_variant))
print_variant_mismatches(var_mismatches_abi3, "Abi3-S209F")

# There is no noticeable difference in gene expression between genotypes so
# we do not validate by expression.


## hAbeta-KI -------------------------------------------------------------------

variants_abeta <- get_variant_mismatches(metadata_all, geno_info,
                                         genotype_pattern = "hAbeta-KI_LoxP_homozygous",
                                         mutation_match = "App-KI",
                                         total_positions = 3)

# Same concept as 5XFAD validation. There are currently no "ambiguous" samples.
variants_merge_abeta <- variants_abeta %>%
  mutate(valid_abeta_ki_variant = (is_carrier & est_genotype == "carrier")  |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_abeta_ki_variant)

valid_samples <- merge(valid_samples, variants_merge_abeta, all = TRUE)
valid_samples$valid_abeta_ki_variant[is.na(valid_samples$valid_abeta_ki_variant)] <- TRUE

var_mismatches_abeta <- merge(variants_abeta, subset(valid_samples, !valid_abeta_ki_variant))
print_variant_mismatches(var_mismatches_abeta, "hAbeta-KI")

# There is no noticeable difference in App expression between genotypes, so we
# do not validate by expression.


## Picalm-H458R ----------------------------------------------------------------

# Variant for Picalm-H458R on the UCI_PrimaryScreen study seems to be
# unreliable, as the variant is being detected in only one of the Picalm samples
# but deletions at that location are found for multiple samples from unrelated
# studies. Therefore we will not validate based on variant calling. It may be
# difficult to detect that variant because it is 3 consecutive base changes. We
# will re-evaluate if a larger Picalm study gets added.

# Additionally, there is a small difference between genotypes (~2 CPM) but not
# enough to be able to confidently validate by expression, so we currently have
# no way to validate the genotype of these mice.

valid_samples$valid_Picalm <- TRUE #"Unknown"


## Trem2-KO --------------------------------------------------------------------

# Trem2-KO isn't reliably detectable from variant calling due to the length of
# base pair deletion in the mutant mice, so we do not validate based on variant
# calling.

# Trem2-KO mice should not express Trem2. Trem2-R47H samples are removed from
# this comparison so only KO and WT samples are compared.
studies_trem2ko <- subset(metadata_all, grepl("Trem2-KO", genotype))
meta_trem2ko <- subset(metadata_all, study_name %in% studies_trem2ko$study_name &
                         !grepl("R47H", genotype))

df_trem2ko <- make_counts_df(meta_trem2ko, counts, symbol_map, "Trem2") %>%
  # Threshold chosen by examination of plot of expression vs genotype
  mutate(expr_trem2 = Trem2 > 7)

valid_trem2ko <- subset(df_trem2ko, (genotype == "Trem2-KO" & !expr_trem2) |
                     (genotype != "Trem2-KO" & expr_trem2))
mismatch_trem2ko <- subset(df_trem2ko, !(merged_file_specimenID %in% valid_trem2ko$merged_file_specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, expr_trem2, Trem2)

print_expression_mismatches(mismatch_trem2ko, "Trem2-KO")

valid_samples$valid_trem2_ko_expression <- !(valid_samples$merged_file_specimenID %in%
                                               mismatch_trem2ko$merged_file_specimenID)


## Trem2-R47H ------------------------------------------------------------------

variants_trem2 <- get_variant_mismatches(metadata_all, geno_info,
                                         genotype_pattern = "Trem2-R47H(_NSS|_CSS)?_(homo|hetero)",
                                         mutation_match = "Trem2-R47H",
                                         total_positions = 1)


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
# carrier mismatches as true mismatches.
variants_merge_trem2 <- variants_trem2 %>%
  mutate(valid_trem2_r47h_variant = is_carrier |
           (!is_carrier & est_genotype != "carrier")) %>%
  select(merged_file_specimenID, valid_trem2_r47h_variant)

valid_samples <- merge(valid_samples, variants_merge_trem2, all = TRUE)
valid_samples$valid_trem2_r47h_variant[is.na(valid_samples$valid_trem2_r47h_variant)] <- TRUE

var_mismatches_trem2 <- merge(variants_trem2, subset(valid_samples, !valid_trem2_r47h_variant))
print_variant_mismatches(var_mismatches_trem2, "Trem2-R47H")

# There is no clear difference in expression of Trem2 between genotypes in each
# study so we can not validate based on expression.


## Bin1 skipped ----------------------------------------------------------------
# No Bin1 carriers present in data

valid_samples$valid_Bin1 <- TRUE #"Unknown"


## Spi1 skipped ----------------------------------------------------------------
# Mutation is present in non-coding part of gene and can't be detected by
# variant calling. There is no discernible difference in expression between
# genotypes so we do not validate by expression.

valid_samples$valid_Spi1 <- TRUE #"Unknown"


# Combine all validation results -----------------------------------------------

valid_columns <- grep("valid_", colnames(valid_samples), value = TRUE)
valid_samples$validated <- rowSums(valid_samples[, valid_columns]) == length(valid_columns)

print(paste(sum(!valid_samples$validated), "samples failed validation:"))
print(subset(valid_samples, !validated))

write.csv(valid_samples, file.path("data", "Model_AD_valid_samples.csv"),
          row.names = FALSE)

# TODO double-check Jax ages match this code
meta_stats <- merge(metadata_all, valid_samples) %>%
  # Fix a few age timepoints for Jax studies
  mutate(
    ageDeath = round(ageDeath),
    ageDeath = case_when(
      study_name == "Jax.IU.Pitt_5XFAD" & ageDeath == 11 ~ 12,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath <= 4 ~ 4,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 4 & ageDeath <= 10 ~ 8,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 10 & ageDeath < 20 ~ 12,
      study_name == "Jax.IU.Pitt_APOE4.Trem2.R47H" & ageDeath > 20 ~ 24,
      .default = ageDeath
    )
  ) %>%
  dplyr::rename(age = ageDeath)

stats <- meta_stats %>%
  group_by(study_name, sex, age, genotype) %>%
  summarize(total_samples = n(),
            valid_samples = sum(validated),
            mismatched_samples = sum(!validated),
            .groups = "drop")

write.csv(stats, "study_validated_samples_stats.csv", row.names = FALSE)

# TODO put on Synapse
