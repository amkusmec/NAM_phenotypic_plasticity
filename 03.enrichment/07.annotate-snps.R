setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(GenomicRanges)


# Data loading and clean-up -----------------------------------------------
### SNPs
snps <- read_delim("~/gxe-gwas/data/geno_2M.map", col_names = FALSE, 
                   delim = "\t", progress = FALSE) %>%
  select(X2, X1, X4)
names(snps) <- c("snpid", "chr", "pos")

### Annotations
source("~/utils/R/get_fgs_intron_exon.R")
fgs <- get_fgs_intron_exon()
genes <- as_data_frame(fgs[["gene"]])
exons <- as_data_frame(fgs[["exon"]])
introns <- as_data_frame(fgs[["intron"]])
rm(fgs)

# Genes
genes <- genes %>%
  mutate(seqname = as.numeric(as.character(seqname))) %>%
  filter(!is.na(seqname)) %>%
  select(seqname, start, end, geneid)
names(genes) <- c("chr", "start", "end", "geneid")

# Exons
exons <- exons %>%
  mutate(seqname = as.numeric(as.character(seqname))) %>%
  filter(!is.na(seqname)) %>%
  mutate(txid = unlist(str_split(txid, "_"))[seq(1, 2*n(), 2)]) %>%
  select(seqname, start, end, txid)
names(exons) <- c("chr", "start", "end", "geneid")

# Introns
introns <- introns %>%
  mutate(seqname = as.numeric(as.character(seqname))) %>%
  filter(!is.na(seqname)) %>%
  mutate(txid = unlist(str_split(txid, "_"))[seq(1, 2*n(), 2)]) %>%
  select(seqname, start, end, txid)
names(introns) <- c("chr", "start", "end", "geneid")

# 5kb upstream of genes
upstream <- genes %>%
  mutate(end = start, start = start - 5000)

# 5kb downstream of genes
downstream <- genes %>%
  mutate(start = end, end = end + 5000)

# HS sites
hs <- bind_rows(read_delim("~/gxe-gwas/data/hs-sites/MNaseHS.Ranges.CrossMap.bed", delim = "\t", col_names = FALSE),
                read_delim("~/gxe-gwas/data/hs-sites/MNaseHS.Root.bed", delim = "\t", col_names = FALSE)) %>%
  mutate(X1 = as.numeric(as.character(X1))) %>%
  filter(!is.na(X1))
names(hs) <- c("chr", "start", "end")


# Create GenomicRanges objects --------------------------------------------
snp.ranges <- with(snps, GRanges(seqnames = chr,
                                 ranges = IRanges(start = pos, end = pos),
                                 snpid = snpid))
exon.ranges <- with(exons, GRanges(seqnames = chr,
                                   ranges = IRanges(start = start, end = end),
                                   geneid = geneid))
intron.ranges <- with(introns, GRanges(seqnames = chr,
                                       ranges = IRanges(start = start, end = end),
                                       geneid = geneid))
upstream.ranges <- with(upstream, GRanges(seqnames = chr,
                                          ranges = IRanges(start = start, end = end),
                                          geneid = geneid))
downstream.ranges <- with(downstream, GRanges(seqnames = chr,
                                              ranges = IRanges(start = start, end = end),
                                              geneid = geneid))
hs.ranges <- with(hs, GRanges(seqnames = chr,
                              ranges = IRanges(start = start, end = end)))
hs.ranges <- reduce(hs.ranges, ignore.strand = TRUE)


# Get the overlaps --------------------------------------------------------
se.hits <- findOverlaps(snp.ranges, exon.ranges, ignore.strand = TRUE)
si.hits <- findOverlaps(snp.ranges, intron.ranges, ignore.strand = TRUE)
su.hits <- findOverlaps(snp.ranges, upstream.ranges, ignore.strand = TRUE)
sd.hits <- findOverlaps(snp.ranges, downstream.ranges, ignore.strand = TRUE)
sh.hits <- findOverlaps(snp.ranges, hs.ranges, ignore.strand = TRUE)


# Combine all the hits in one data frame ----------------------------------
all.indices <- c(su.hits@from, se.hits@from, si.hits@from, sd.hits@from, sh.hits@from)
anno <- data_frame(snpid = snp.ranges@elementMetadata@listData$snpid[all.indices],
                   type = c(rep("5kb upstream", length(su.hits@from)),
                            rep("Exon", length(se.hits@from)),
                            rep("Intron", length(si.hits@from)),
                            rep("5kb downstream", length(sd.hits@from)),
                            rep("Hypersensitivity site", length(sh.hits@from))),
                   geneid = c(upstream.ranges@elementMetadata@listData$geneid[su.hits@to],
                              exon.ranges@elementMetadata@listData$geneid[se.hits@to],
                              intron.ranges@elementMetadata@listData$geneid[si.hits@to],
                              downstream.ranges@elementMetadata@listData$geneid[sd.hits@to],
                              rep("NA", length(sh.hits@from))))

anno <- full_join(anno, snps, by = "snpid")
anno$type[is.na(anno$type)] <- "Intergenic"

anno <- anno %>%
  mutate(key = paste0(snpid, type, geneid)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(snpid, chr:pos, type:geneid)


# Set-up the type factor --------------------------------------------------
# Modified from Rodgers-Melnick et al. (2016)
anno$type <- factor(anno$type, ordered = TRUE,
                    levels = c("Exon", "5kb upstream", "5kb downstream",
                               "Intron", "Hypersensitivity site", "Intergenic"))
anno <- anno %>%
  arrange(chr, pos, type) %>%
  group_by(snpid) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Save the result
saveRDS(anno, "data/snp-annotations.rds")
