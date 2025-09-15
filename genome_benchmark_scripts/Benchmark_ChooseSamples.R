# This script picks samples that will be used to benchmark custom reference
# genomes. Samples are selected from three studies: Jax.IU.Pitt_5XFAD,
# Jax.IU.Pitt_APOE4.Trem2.R47H, UCI_3xTg-AD, and UCI_Clu-h2kbKI, which cover all
# human transgenes present in any MODEL-AD study. Samples were selected from
# female mice that were approximately 12 months of age, allowing for differences
# in sample availability:
#   3xTg: 12 months
#   5XFAD: 11 months
#   APOE4: 13-16 months
#   Clu: 12 months
#
# The scripts 01_Harmonize_Metadata.R and 02_NF_CreateSampleSheets.R must be run
# before running this script.

library(synapser)
library(stringr)
library(dplyr)

synLogin()
tmp_dir <- file.path("data", "tmp")
dir.create(tmp_dir, showWarnings = FALSE)

meta_file <- synGet("syn61850266", downloadLocation = tmp_dir)
metadata <- read.csv(file.path(tmp_dir, "Model_AD_merged_metadata.csv"))

# Remove samples that have an NA genotype
metadata <- subset(metadata, !is.na(genotype))

# 5XFAD, APOE4, 3xTg, and Clu cover all the humanized genes in the custom reference genome
bench <- subset(metadata, study_name %in% c("Jax.IU.Pitt_5XFAD",
                                            "Jax.IU.Pitt_APOE4.Trem2.R47H",
                                            "UCI_3xTg-AD",
                                            "UCI_Clu-h2kbKI"))

# Select female mice, and narrow down the ages to the ranges we want
bench <- subset(bench, sex == "female" &
                  round(ageDeath) %in% c(11, 12, 13, 14, 15, 16))

# Get rid of some genotypes that only have 1-2 samples
low_genotypes <- c("APOE4-KI_heterozygous; Trem2-R47H_homozygous",
                   "APOE4-KI_homozygous; Trem2-R47H_heterozygous")
bench <- subset(bench, !(genotype %in% low_genotypes))

# Check -- all specimen IDs should be unique
stopifnot(length(unique(bench$specimenID)) == length(bench$specimenID))

# Write this subset of metadata to a csv file in case we need it later
write.csv(bench,
          file.path("data", "genome_reference_benchmark_metadata.csv"),
          row.names = FALSE, quote = FALSE)

# Assuming we have already run 02_NF_CreateSampleSheets for all of these
# studies, we can load them in and filter to just the benchmark samples
lines <- lapply(unique(bench$study_name), function(study) {
  ss_name <- file.path("data", "sample_sheets",
                       paste0(study, "_samplesheet_rnaseq.csv"))
  stopifnot(file.exists(ss_name))

  samples <- read.csv(ss_name)
  samples <- subset(samples, sample %in% bench$specimenID)
  return(samples)
})

lines <- do.call(rbind, lines)

write.csv(lines,
          file.path("data", "sample_sheets", "genome_reference_samplesheet_rnaseq.csv"),
          row.names = FALSE, quote = FALSE)
