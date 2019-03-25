source("00.load-packages.R")
library(GenomicRanges)


# Prep the significant SNPs -----------------------------------------------
sig <- readRDS("gwas-results-iqr/sig-snps.rds")
mranges <- sig %>%
  filter(Type == "MainEffect") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position, end = Position),
                  pheno = Phenotype, type = Type, SNP = SNP))
sranges <- sig %>%
  filter(Type == "Slope") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position, end = Position),
                  pheno = Phenotype, type = Type, SNP = SNP))
vranges <- sig %>%
  filter(Type == "VarE") %>%
  with(., GRanges(seqnames = Chromosome,
                  ranges = IRanges(start = Position, end = Position),
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


# Control and storage variables -------------------------------------------
windows <- c(2e3, 5e3, 1e4, 1.5e4, 2e4, 3e4, 4e4, 5e4)
overlaps <- vector(mode = "numeric", length = length(windows))
totals <- vector(mode = "numeric", length = length(windows))


# Get overlaps at varying window sizes ------------------------------------
for (i in seq_along(windows)) {
  # Resize the SNP windows
  mranges <- resize(mranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  sranges <- resize(sranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  vranges <- resize(vranges, width = 2*windows[i] + 1, fix = "center", use.names = FALSE)
  
  # Get the candidate genes
  mgenes <- findOverlaps(geneRanges, mranges, ignore.strand = TRUE)
  sgenes <- findOverlaps(geneRanges, sranges, ignore.strand = TRUE)
  vgenes <- findOverlaps(geneRanges, vranges, ignore.strand = TRUE)
  
  # Index variables for concision
  gmind <- mgenes@from; mind <- mgenes@to
  gsind <- sgenes@from; sind <- sgenes@to
  gvind <- vgenes@from; vind <- vgenes@to
  
  # Create a results table
  d <- data_frame(GeneID = geneRanges@elementMetadata@listData$ids[c(gmind, gsind, gvind)],
                  SNP = c(mranges@elementMetadata@listData$SNP[mind],
                          sranges@elementMetadata@listData$SNP[sind],
                          vranges@elementMetadata@listData$SNP[vind]),
                  Phenotype = c(mranges@elementMetadata@listData$pheno[mind],
                                sranges@elementMetadata@listData$pheno[sind],
                                vranges@elementMetadata@listData$pheno[vind]),
                  Type = c(mranges@elementMetadata@listData$type[mind],
                           sranges@elementMetadata@listData$type[sind],
                           vranges@elementMetadata@listData$type[vind]))
  
  # Save the results
  totals[i] <- length(unique(d$GeneID))
  overlaps[i] <- d %>% count(Phenotype, GeneID) %>% filter(n > 1) %>% nrow(.)
}


# Summarize the results ---------------------------------------------------
d <- data_frame(WindowSize = rep(2*windows/1000, times = 2),
                Genes = c(totals, overlaps),
                Facet = rep(c("Total", "Overlap"), each = length(windows)))
saveRDS(d, "gwas-results-iqr/windows.rds")

f4b <- ggplot(d, aes(x = WindowSize, y = Genes)) + geom_point() + geom_line() + theme_bw() +
  facet_wrap("Facet", nrow = 2, scales = "free_y", strip.position = "left") +
  labs(y = "Number of Candidate Genes", x = "Window Size (kb)") +
  scale_x_continuous(breaks = 2*windows/1000, labels = 2*windows/1000) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
# ggsave("figures/windows.png", width = 5, height = 4, units = "in", dpi = 300)
