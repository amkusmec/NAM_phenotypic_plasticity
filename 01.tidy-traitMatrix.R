setwd("~/gxe-gwas2")

library(tibble)
library(readr)
library(tidyr)
library(dplyr)
library(magrittr)
library(stringr)

# Load the rownames (genotypes)
# The first line is a header, so remove it
taxa <- readLines("data/nam282_traitMatrix_rownames.txt")[-1]

# Load the full trait matrix
traitMatrix <- read_csv("data/nam282_traitMatrix_filtered.txt", 
                        col_names = TRUE, na = "NaN") %>%
  mutate(Genotype = taxa)

# Full trait matrix is 5,706 taxa by 254 phenotype-environment combinations
# Tidy the data by splitting the phenotype-environment combinations and putting
# all measurements in a Measure column
traitMatrix <- traitMatrix %>%
  gather(key = Phenotype, value = Measure, CobDiameter_id065:VarName_id08A) %>%
  separate(Phenotype, c("Phenotype", "Environment"), sep = "_id") %>%
  select(Genotype, Phenotype, Environment, Measure)

# Table has 1,443,618 observations at this point
# Keep only non-missing observations from US-NAM for genotypes for which
# we have SNPs and only certain phenotypes
genos <- readLines("data/genotypes_to_keep.txt")
phenos <- readLines("data/phenotypes_to_keep.txt")
traitMatrix <- traitMatrix %>%
  filter(Genotype %in% genos, Phenotype %in% phenos, !is.na(Measure))

# Construct a genotype-phenotype key
# Remove all keys that do not occur at least 3 times
traitMatrix <- traitMatrix %>% mutate(Key = paste0(Genotype, Phenotype))
keyCounts <- traitMatrix %>%
  dplyr::count(Key) %>%
  filter(n >= 3)
traitMatrix <- traitMatrix %>%
  filter(Key %in% keyCounts$Key) %>%
  select(-Key)

# Final dataset contains 749,590 observations
saveRDS(traitMatrix, "data/tidy_traitMatrix.rds")
