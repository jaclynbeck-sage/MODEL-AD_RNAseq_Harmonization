library(synapser)
library(dplyr)

# Set up -----------------------------------------------------------------------

# Synapse IDs used in this script -- merged metadata and the transcript_counts
# file with data from all studies
syn_metadata_file_id <- "syn61850266.3"
syn_gene_counts_id <- "syn62690577.1"
syn_symbol_map_id <- "syn62063692.1"

synLogin()

github <- "https://github.com/jaclynbeck-sage/MODEL-AD_RNAseq_Harmonization/blob/main/05_Sample_Validation.Rmd"
tmp_dir <- file.path("data", "tmp")


# Load counts and metadata -----------------------------------------------------

metadata_file <- synGet(syn_metadata_file_id, downloadLocation = tmp_dir)
counts_file <- synGet(syn_gene_counts_id, downloadLocation = tmp_dir)
symbol_map_file <- synGet(syn_symbol_map_id, downloadLocation = tmp_dir)

metadata_all <- read.csv(metadata_file$path)
counts <- read.table(counts_file$path, header = TRUE, sep = "\t",
                     row.names = 1) %>%
  # Get rid of transcript_id column
  select(-transcript_id.s.)

# Not all samples in the metadata file appear in the counts matrix and vice versa
metadata_all <- subset(metadata_all, merged_file_specimenID %in% colnames(counts)) %>%
  select(individualID, specimenID, merged_file_specimenID, sex, genotype, study_name)
counts <- counts[, metadata_all$merged_file_specimenID]

symbol_map <- read.csv(symbol_map_file$path)

# Convert counts to CPM so samples are comparable to each other
counts <- sweep(counts, 2, colSums(counts), "/") * 1e6


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

# This will be combined with the genotype-validated samples
sex_matched_samples <- subset(xy_df, valid_sex == TRUE)
sex_matched_samples <- sex_matched_samples$merged_file_specimenID

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

genotype_matched_samples <- c()
mismatched_samples <- c()

## 5XFAD carrier / non-carriers ------------------------------------------------

meta_5x <- subset(metadata_all, grepl("5XFAD", genotype))

df_5x <- make_counts_df(meta_5x, counts, symbol_map,
                        c("APP", "PSEN1")) %>%
  # Using > 1 CPM is sufficient for this genotype
  mutate(est_5x = APP > 1 & PSEN1 > 1)

valid_5x <- subset(df_5x, (grepl("5XFAD_carrier", genotype) & est_5x) |
                     (grepl("5XFAD_noncarrier", genotype) & !est_5x))
mismatch_5x <- subset(df_5x, !(specimenID %in% valid_5x$specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, est_5x, APP, PSEN1)

if (nrow(mismatch_5x) > 0) {
  print(paste(nrow(mismatch_5x), "samples have a 5XFAD genotype mismatch:"))
  print(select(mismatch_5x, -merged_file_specimenID))
} else {
  print("No samples have a 5XFAD genotype mismatch.")
}

genotype_matched_samples <- c(genotype_matched_samples, valid_5x$merged_file_specimenID)
mismatched_samples <- c(mismatched_samples, mismatch_5x$merged_file_specimenID)


## 3xTg carrier / non-carriers -------------------------------------------------

meta_3x <- subset(metadata_all, grepl("3xTg-AD", genotype))

df_3x <- make_counts_df(meta_3x, counts, symbol_map,
                        c("APP", "MAPT")) %>%
  # Using > 1 CPM is sufficient for this genotype
  mutate(est_3x = APP > 1 & MAPT > 1)

valid_3x <- subset(df_3x, (grepl("3xTg-AD_carrier", genotype) & est_3x) |
                     (grepl("3xTg-AD_noncarrier", genotype) & !est_3x))
mismatch_3x <- subset(df_3x, !(specimenID %in% valid_3x$specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, est_3x, APP, MAPT)

if (nrow(mismatch_3x) > 0) {
  print(paste(nrow(mismatch_3x), "samples have a 3xTg-AD genotype mismatch:"))
  print(select(mismatch_3x, -merged_file_specimenID))
} else {
  print("No samples have a 3xTg-AD genotype mismatch.")
}

genotype_matched_samples <- c(genotype_matched_samples, valid_3x$merged_file_specimenID)
mismatched_samples <- c(mismatched_samples, mismatch_3x$merged_file_specimenID)


## APOE4-KI homozygous / noncarrier --------------------------------------------

meta_apoe <- subset(metadata_all, grepl("APOE4-KI", genotype))

df_apoe <- make_counts_df(meta_apoe, counts, symbol_map,
                        c("APOE")) %>%
  # Using > 1 CPM is sufficient for this genotype
  mutate(est_apoe = APOE > 1)

valid_apoe <- subset(df_apoe, (grepl("APOE4-KI_(homo|hetero)zygous", genotype) & est_apoe) |
                     (grepl("APOE4-KI_noncarrier", genotype) & !est_apoe))
mismatch_apoe <- subset(df_apoe, !(specimenID %in% valid_apoe$specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, est_apoe, APOE)

if (nrow(mismatch_apoe) > 0) {
  print(paste(nrow(mismatch_apoe), "samples have an APOE4-KI genotype mismatch:"))
  print(select(mismatch_apoe, -merged_file_specimenID))
} else {
  print("No samples have an APOE4-KI genotype mismatch.")
}

genotype_matched_samples <- c(genotype_matched_samples, valid_apoe$merged_file_specimenID)
mismatched_samples <- c(mismatched_samples, mismatch_apoe$merged_file_specimenID)

## Everything else -------------------------------------------------------------
# None of these should express any human genes

no_hu_meta <- subset(metadata_all,
                     !(merged_file_specimenID %in% genotype_matched_samples) &
                       !(merged_file_specimenID %in% mismatched_samples))

no_hu_df <- make_counts_df(no_hu_meta, counts, symbol_map,
                           c("APOE", "APP", "MAPT", "PSEN1")) %>%
  # Some low-level counts are possible but > 1 CPM shouldn't happen for
  # genotypes that don't contain transgenes
  mutate(hu_expr = APOE > 1 | APP > 1 | MAPT > 1 | PSEN1 > 1)

valid_no_hu <- subset(no_hu_df, !hu_expr)
mismatch_no_hu <- subset(no_hu_df, !(specimenID %in% valid_no_hu$specimenID)) %>%
  select(study_name, specimenID, merged_file_specimenID, genotype, hu_expr, APOE:PSEN1)

if (nrow(mismatch_no_hu) > 0) {
  print(paste(nrow(mismatch_no_hu), "samples erroneously express transgenes:"))
  print(select(mismatch_no_hu, -merged_file_specimenID))
} else {
  print("No samples erroneously express transgenes.")
}

genotype_matched_samples <- c(genotype_matched_samples, valid_no_hu$merged_file_specimenID)
mismatched_samples <- c(mismatched_samples, mismatch_no_hu$merged_file_specimenID)

# Make sure no genotype matched samples ended up being mismatches in a different
# comparison
genotype_matched_samples <- setdiff(genotype_matched_samples, mismatched_samples)


# Combine sex and genotype results ---------------------------------------------

meta_validated <- metadata_all %>%
  mutate(valid_sex = merged_file_specimenID %in% sex_matched_samples,
         valid_genotype = merged_file_specimenID %in% genotype_matched_samples,
         is_valid = valid_sex & valid_genotype)

print("Valid samples per study:")
print(table(meta_validated$study_name, meta_validated$is_valid))

write.csv(meta_validated, file.path("data", "Model_AD_merged_metadata_valid_samples.csv"))

# TODO put on Synapse




