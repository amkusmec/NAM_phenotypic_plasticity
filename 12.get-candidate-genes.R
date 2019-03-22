setwd("~/gxe-gwas2/gwas-results-iqr")

source("../00.load-packages.R")
library(GenomicRanges)


# Prep the significant SNPs -----------------------------------------------
sig <- readRDS("sig-snps.rds")
window <- 1e4

mranges <- sig %>%
  filter(Type == "MainEffect") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position - window, end = Position + window),
                  pheno = Phenotype, type = Type, SNP = SNP))
sranges <- sig %>%
  filter(Type == "Slope") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position - window, end = Position + window),
                  pheno = Phenotype, type = Type, SNP = SNP))
vranges <- sig %>%
  filter(Type == "VarE") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position - window, end = Position + window),
                  pheno = Phenotype, type = Type, SNP = SNP))


# Prep the genes ----------------------------------------------------------
gff <- read_delim("~/anno/ZmB73_5b_FGS.gff", comment = "#", delim = "\t", 
                  progress = FALSE, col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  select(X1, X4, X5, X9) %>%
  mutate(X9 = unlist(str_split(X9, ";"))[seq(1, 3*n(), 3)]) %>%
  mutate(X9 = gsub('ID=', '', X9)) %>%
  arrange(X1, X4)
names(gff) <- c("chr", "start", "end", "gene")

geneRanges <- with(gff, GRanges(seqnames = chr,
                                ranges = IRanges(start = start, end = end),
                                ids = gene))


# Get the candidate genes -------------------------------------------------
mgenes <- findOverlaps(geneRanges, mranges, ignore.strand = TRUE)
sgenes <- findOverlaps(geneRanges, sranges, ignore.strand = TRUE)
vgenes <- findOverlaps(geneRanges, vranges, ignore.strand = TRUE)


# Create a results table --------------------------------------------------
# Index variables for concision
gmind <- mgenes@from; mind <- mgenes@to
gsind <- sgenes@from; sind <- sgenes@to
gvind <- vgenes@from; vind <- vgenes@to

d <- data_frame(GeneID = geneRanges@elementMetadata@listData$ids[c(gmind, gsind, gvind)],
                SNP = c(mranges@elementMetadata@listData$SNP[mind],
                        sranges@elementMetadata@listData$SNP[sind],
                        vranges@elementMetadata@listData$SNP[vind]),
                Phenotype = c(mranges@elementMetadata@listData$pheno[mind],
                              sranges@elementMetadata@listData$pheno[sind],
                              vranges@elementMetadata@listData$pheno[vind]),
                Type = c(mranges@elementMetadata@listData$type[mind],
                         sranges@elementMetadata@listData$type[sind],
                         vranges@elementMetadata@listData$type[vind])) %>%
  inner_join(., gff, by = c("GeneID" = "gene")) %>%
  select(SNP:Type, GeneID, chr:end) %>%
  arrange(Phenotype, Type, chr, start)

# Save a version for downstream analyses
saveRDS(d, "candidate-genes.rds")

# Save a version for supplementary tables
names(d) <- c("SNP", "Phenotype", "Measure", "Gene", "Chromosome", "Start", "End")
d %>% mutate(Measure = gsub("MainEffect", "Mean phenotype", Measure),
             Measure = gsub("Slope", "Linear plasticity", Measure),
             Measure = gsub("VarE", "Non-linear plasticity", Measure)) %>%
  write_csv(., path = "../supplementary/candidate_genes.csv")
