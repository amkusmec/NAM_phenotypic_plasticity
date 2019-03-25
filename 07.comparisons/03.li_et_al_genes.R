setwd("~/gxe-gwas2")

library(GenomicRanges)
cv <- read.csv("~/gxe-gwas/data/floweringcv_qtl.csv")
cv <- subset(cv, DACV == "_" | DSCV == "_")
cv$upper <- 1e6*cv$upper
cv$lower <- 1e6*cv$lower
# upper and lower are swapped in the file
cv_ranges <- with(cv, GRanges(seqnames = chr,
                              ranges = IRanges(start = upper, end = lower)))

gff <- read.delim("~/anno/ZmB73_5b_FGS.gff", comment.char = "#", header = FALSE)
gff <- gff[gff$V3 == "gene", c(1, 4, 5, 9)]
gff$V9 <- as.character(gff$V9)
gff$V9 <- unlist(str_split(gff$V9, ";"))[seq(1, 3*nrow(gff), 3)]
gff$V9 <- gsub('ID=', '', gff$V9)
gff$V1 <- as.character(gff$V1)
gff$V1 <- as.numeric(gff$V1)
gff <- gff[!is.na(gff$V1), ]
names(gff) <- c("chr", "start", "end", "gene")
gff <- gff[order(gff$chr, gff$start), ]
geneRanges <- with(gff, GRanges(seqnames = chr,
                                ranges = IRanges(start = start, end = end),
                                ids = gene))

cv_overlaps <- findOverlaps(geneRanges, cv_ranges, ignore.strand = TRUE)
length(unique(cv_overlaps@from)) # 6,197 genes in these intervals

cv$pos <- 1e6*cv$pos
cv_ranges_small <- with(cv, GRanges(seqnames = chr,
                                    ranges = IRanges(start = pos - 1e6, end = pos + 1e6)))
cv_overlaps_small <- findOverlaps(geneRanges, cv_ranges_small, ignore.strand = TRUE)
length(unique(cv_overlaps_small@from)) # 1,848 genes in 1 Mb intervals centered on peaks

cv_overlaps_counts <- sapply(split(cv_overlaps@from, cv_overlaps@to), length, simplify = TRUE)
summary(cv_overlaps_counts)

cv_overlaps_small_counts <- sapply(split(cv_overlaps_small@from, cv_overlaps_small@to), length, simplify = TRUE)
summary(cv_overlaps_small_counts)
