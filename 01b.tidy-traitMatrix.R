setwd("~/gxe-gwas2")

source("00.load-packages.R")


# Load the trait matrix ---------------------------------------------------
traitMatrix <- read_tsv("data/traitMatrix_maize282NAM_v15-130212.txt",
                        col_names = FALSE, skip = 2)
phenos <- read_lines("data/traitMatrix_maize282NAM_v15-130212.txt", n_max = 2)
envs <- unlist(str_split(phenos[2], "\t"))[-1]
phenos <- unlist(str_split(phenos[1], "\t"))[-1]
phenos <- gsub("-", "", phenos)
phenos <- paste(phenos, envs, sep = "_")
phenos <- c("Taxa", phenos)
names(traitMatrix) <- phenos


# Tidy the data -----------------------------------------------------------
# Full trait matrix is 5,706 taxa by 254 phenotype-environment combinations
# Tidy the data by splitting the phenotype-environment combinations and putting
# all measurements in a Measure column
traitMatrix <- traitMatrix %>%
  gather(key = Phenotype, value = Measure, 
         GerminationCount_065:HerbicideSensitivity_08A) %>%
  separate(Phenotype, c("Phenotype", "Environment"), sep = "_") %>%
  rename(Genotype = Taxa) %>%
  select(Genotype, Phenotype, Environment, Measure)


# Further filtering -------------------------------------------------------
# Table has 1,626,210 observations at this point
# Keep only non-missing observations from US-NAM for genotypes for which
# we have SNPs and only certain phenotypes
genos <- readLines("data/genotypes_to_keep.txt")
phenos <- readLines("data/phenotypes_to_keep.txt")
traitMatrix <- traitMatrix %>%
  filter(Genotype %in% genos, Phenotype %in% phenos, !is.na(Measure))


# IQR Filtering -----------------------------------------------------------
# 754,368 observations
# Construct cutoffs based on +/- 1.5*IQR RIL- and Environment-wise
ril_cutoff <- traitMatrix %>%
  group_by(Phenotype, Genotype) %>%
  summarise(R_Q1 = quantile(Measure, probs = 0.25),
            R_Q3 = quantile(Measure, probs = 0.75)) %>%
  ungroup() %>%
  mutate(R_IQR = R_Q3 - R_Q1,
         R_Lower = R_Q1 - 1.5*R_IQR,
         R_Upper = R_Q3 + 1.5*R_IQR)
env_cutoff <- traitMatrix %>%
  group_by(Phenotype, Environment) %>%
  summarise(E_Q1 = quantile(Measure, probs = 0.25),
            E_Q3 = quantile(Measure, probs = 0.75)) %>%
  ungroup() %>%
  mutate(E_IQR = E_Q3 - E_Q1,
         E_Lower = E_Q1 - 1.5*E_IQR,
         E_Upper = E_Q3 + 1.5*E_IQR)
traitMatrix <- traitMatrix %>%
  inner_join(., ril_cutoff, by = c("Phenotype", "Genotype")) %>%
  inner_join(., env_cutoff, by = c("Phenotype", "Environment"))
traitMatrix <- traitMatrix %>%
  filter(Measure >= R_Lower, Measure <= R_Upper,
         Measure >= E_Lower, Measure <= E_Upper) %>%
  select(Genotype:Measure)


# Construct a genotype-phenotype key --------------------------------------
# 701,728 observations
# Remove all keys that do not occur at least 3 times
traitMatrix <- traitMatrix %>% mutate(Key = paste0(Genotype, Phenotype))
keyCounts <- traitMatrix %>%
  count(Key) %>%
  filter(n >= 3)
traitMatrix <- traitMatrix %>%
  filter(Key %in% keyCounts$Key) %>%
  select(-Key)

# Final dataset contains 698,738 observations
saveRDS(traitMatrix, "data/tidy_traitMatrix_IQR_AK.rds")
