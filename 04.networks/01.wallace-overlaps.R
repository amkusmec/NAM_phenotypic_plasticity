setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(GenomicRanges)
library(scales)


# Common phenotypes between the two datasets (14/23) ---------------------------
keep <- c("Boxcox-transformed leaf angle", "Cob diameter", "Days to silk", "Days to anthesis",
          "Ear row number", "Ear height", "Height above ear", "Leaf length", "Leaf width",
          "Photoperiod growing-degree days to anthesis", "Photoperiod Growing-degree days to silk",
          "Plant height", "Tassel branch number", "Tassel length")
mine <- c("UpperLeafAngle", "CobDiameter", "DaystoSilk", "DaysToTassel",
          "EarRowNumber", "EarHeight", "HeightAboveEar", "LeafLength", "LeafWidth",
          "GDDDaystoTassel", "GDDDaystoSilk",
          "PlantHeight", "TasselPrimaryBranches", "TasselLength")


# Prep the Wallace dataset ------------------------------------------------
wallace <- read_delim("~/gxe-gwas/data/Wallace_etal_2014_PLoSGenet_GWAS_hits-150112.txt",
                      delim = "\t") %>%
  filter(trait %in% keep)
for (i in seq_along(keep)) {
  wallace$trait <- gsub(keep[i], mine[i], wallace$trait)
}


# Prep the SNPs -----------------------------------------------------------
sig <- readRDS("gwas-results-iqr/sig-snps.rds") %>%
  filter(Phenotype %in% mine)


# Construct GenomicRanges object ------------------------------------------
wranges <- with(wallace, GRanges(seqnames = chr,
                                 ranges = IRanges(start = pos, end = pos),
                                 pheno = trait,
                                 SNP = paste(chr, pos, sep = "_")))
mranges <- with(subset(sig, Type == "MainEffect"),
                GRanges(seqnames = Chromosome,
                        ranges = IRanges(start = Position, end = Position),
                        pheno = Phenotype,
                        type = Type,
                        SNP = SNP))
sranges <- with(subset(sig, Type == "Slope"),
                GRanges(seqnames = Chromosome,
                        ranges = IRanges(start = Position, end = Position),
                        pheno = Phenotype,
                        type = Type,
                        SNP = SNP))
vranges <- with(subset(sig, Type == "VarE"),
                GRanges(seqnames = Chromosome,
                        ranges = IRanges(start = Position, end = Position),
                        pheno = Phenotype,
                        type = Type,
                        SNP = SNP))


# Window sizes to test ----------------------------------------------------
windows <- c(2e3, 5e3, 1e4, 1.5e4, 2e4, 3e4, 4e4, 5e4)
types <- c("MainEffect", "Slope", "VarE")
overlaps <- data_frame(WindowSize = rep(windows, times = length(types)),
                       Type = rep(types, each = length(windows)),
                       Overlaps = 0)


# Get the number of overlaps ----------------------------------------------
for (i in seq_along(windows)) {
  # Resize the SNP ranges
  mranges <- resize(mranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  sranges <- resize(sranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  vranges <- resize(vranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  wranges <- resize(wranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  
  # Get the overlaps
  msnps <- findOverlaps(wranges, mranges, ignore.strand = TRUE)
  ssnps <- findOverlaps(wranges, sranges, ignore.strand = TRUE)
  vsnps <- findOverlaps(wranges, vranges, ignore.strand = TRUE)
  
  # Indices for subsetting
  wmind <- msnps@from; mind <- msnps@to
  wsind <- ssnps@from; sind <- ssnps@to
  wvind <- vsnps@from; vind <- vsnps@to
  
  # Create a table of overlapping bins with associated phenotypes
  d <- data_frame(MySNP = c(mranges@elementMetadata@listData$SNP[mind],
                            sranges@elementMetadata@listData$SNP[sind],
                            vranges@elementMetadata@listData$SNP[vind]),
                  TheirSNP = wranges@elementMetadata@listData$SNP[c(wmind, wsind, wvind)],
                  MyPheno = c(mranges@elementMetadata@listData$pheno[mind],
                              sranges@elementMetadata@listData$pheno[sind],
                              vranges@elementMetadata@listData$pheno[vind]),
                  TheirPheno = wranges@elementMetadata@listData$pheno[c(wmind, wsind, wvind)],
                  Type = c(mranges@elementMetadata@listData$type[mind],
                           sranges@elementMetadata@listData$type[sind],
                           vranges@elementMetadata@listData$type[vind]),
                  WindowSize = windows[i])
  
  # Count the overlaps
  for (j in seq_along(types)) {
    overlaps[overlaps$WindowSize == windows[i] & overlaps$Type == types[j], "Overlaps"] <-
      with(subset(d, Type == types[j]), sum(MyPheno == TheirPheno))
  }
}


# Load the results of permutation testing ---------------------------------
# Calculate 2.5 and 97.5 percentiles
emp <- bind_rows(as_data_frame(readRDS("data/overlap-permutations-iqr-MainEffect.rds")),
                 as_data_frame(readRDS("data/overlap-permutations-iqr-Slope.rds")),
                 as_data_frame(readRDS("data/overlap-permutations-iqr-VarE.rds"))) %>%
  mutate(Type = rep(types, each = 1000)) %>%
  gather(key = WindowSize, value = Overlaps, X2000:X50000) %>%
  group_by(Type, WindowSize) %>%
  summarise(X2.5 = quantile(Overlaps, probs = 0.025),
            X5 = quantile(Overlaps, probs = 0.05),
            X95 = quantile(Overlaps, probs = 0.95),
            X97.5 = quantile(Overlaps, probs = 0.975)) %>%
  mutate(WindowSize = as.integer(gsub("X", "", WindowSize))) %>%
  arrange(Type, WindowSize)

# Determine which overlaps are significantly different from the null
overlaps <- overlaps %>%
  full_join(., emp, by = c("Type", "WindowSize")) %>%
  mutate(Sig2 = if_else(Overlaps <= X2.5 | Overlaps >= X97.5, TRUE, FALSE),
         SigL = if_else(Overlaps <= X5, TRUE, FALSE),
         SigU = if_else(Overlaps >= X95, TRUE, FALSE))


# Format for plotting -----------------------------------------------------
overlaps <- overlaps %>%
  mutate(WindowSize = 2*WindowSize/1000,
         Type = gsub("MainEffect", "Mean Phenotype", Type),
         Type = gsub("Slope", "Linear Plasticity", Type),
         Type = gsub("VarE", "Non-linear Plasticity", Type),
         Type = factor(Type, ordered = TRUE,
                       levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")))
overlaps$Overlaps[1:8] <- overlaps$Overlaps[1:8]/length(unique(subset(sig, Type == "MainEffect")$SNP))
overlaps$Overlaps[9:16] <- overlaps$Overlaps[9:16]/length(unique(subset(sig, Type == "Slope")$SNP))
overlaps$Overlaps[17:24] <- overlaps$Overlaps[17:24]/length(unique(subset(sig, Type == "VarE")$SNP))


# Plot --------------------------------------------------------------------
f3b <- ggplot(overlaps, aes(x = WindowSize, y = Overlaps, colour = Type)) +
  geom_line(size = 1) + theme_bw() + guides(shape = FALSE) +
  geom_point(data = subset(overlaps, Sig2), aes(colour = Type), shape = 16, size = 3.5) +
  geom_point(data = subset(overlaps, !Sig2), aes(colour = Type), shape = 1, size = 3.5) +
  labs(y = "% of Overlapping Windows", x = "Window Size (kb)", colour = "") +
  scale_x_continuous(breaks = unique(overlaps$WindowSize), 
                     labels = unique(overlaps$WindowSize)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = colors) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  guides(colour = "none")
        # legend.position = "top",
        # legend.text = element_text(size = 8))
# ggsave("figures/wallace-overlaps.png", width = 4.5, height = 4.5, units = "in", dpi = 300)

# saveRDS(overlaps, "gwas-results/wallace-overlaps.rds")
