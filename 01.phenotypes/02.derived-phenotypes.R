### Calculate two derived phenotypes:
###  1) KernelDepth = EarDiameter - CobDiameter
###  2) HeightAboveEar = PlantHeight - EarHeight

source("00.load-packages.R")

# Load the trait matrix
traitMatrix <- readRDS("data/tidy_traitMatrix_IQR_AK.rds")

# Calculate the kernel depth phenotype
ed <- subset(traitMatrix, Phenotype == "EarDiameter")
cd <- subset(traitMatrix, Phenotype == "CobDiameter")
temp <- merge(ed, cd, by = c("Genotype", "Environment"), all = FALSE, sort = FALSE)
temp <- temp %>%
  mutate(Measure = (Measure.x - Measure.y)/2,
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

# Order everything for convenience
traitMatrix <- traitMatrix %>%
  arrange(Phenotype, Environment, Genotype)

# Final dataset has dimensions
#  - 4,890 genotypes
#  - 12 environments
#  - 23 phenotypes
#  - 766,300 total observations
saveRDS(traitMatrix, "data/tidy_traitMatrix_IQR_AK.rds")


# Supplementary trait distributions ---------------------------------------
traitMatrix <- split(traitMatrix, traitMatrix$Phenotype)
traitMatrix %>% map(function(x) {
  p <- ggplot(x, aes(x = Environment, y = Measure)) + geom_boxplot() + 
    theme_bw() + ggtitle(x$Phenotype[1])
  ggsave(paste0("supplementary/dist_", x$Phenotype[1], ".png"), plot = p, 
         height = 5, width = 5, dpi = 200)
  NULL
})
