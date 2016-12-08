setwd("~/gxe-gwas2")

# Load the trait matrix
traitMatrix <- readRDS("data/tidy_traitMatrix.rds")

# Calculate the kernel depth phenotype
ed <- subset(traitMatrix, Phenotype == "EarDiameter")
cd <- subset(traitMatrix, Phenotype == "CobDiameter")
temp <- merge(ed, cd, by = c("Genotype", "Environment"), all = FALSE, sort = FALSE)
temp <- temp %>%
  mutate(Measure = Measure.x - Measure.y,
         Phenotype = "KernelDepth") %>%
  select(Genotype, Phenotype, Environment, Measure)
traitMatrix <- rbind(traitMatrix, temp)

# Calculate the height above ear phenotype
ph <- subset(traitMatrix, Phenotype == "PlantHeight")
eh <- subset(traitMatrix, Phenotype == "EarHeight")
temp <- merge(ph, eh, by = c("Genotype", "Environment"), all = FALSE, sort = FALSE)
temp <- temp %>%
  mutate(Measure = Measure.x - Measure.y,
         Phenotype = "HeightAboveEar") %>%
  select(Genotype, Phenotype, Environment, Measure)
traitMatrix <- rbind(traitMatrix, temp)

# Final dataset has dimensions
#  - 4,890 genotypes
#  - 12 environments
#  - 23 phenotypes
#  - 825,969 total observations
saveRDS(traitMatrix, "data/tidy_traitMatrix.rds")
