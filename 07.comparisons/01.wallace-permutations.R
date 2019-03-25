setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(GenomicRanges)
library(doParallel)
library(foreach)
library(iterators)
library(argparse)


# Set-up the parser -------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("-t", "--type", type = "character")
args <- parser$parse_args()
t <- args$type


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
snps <- readRDS("data/snp-annotations.rds")
snps$snpid <- paste0("X", snps$snpid)
names(snps)[4] <- "Location"
sig <- readRDS("gwas-results-iqr/sig-snps.rds")
sig <- left_join(sig, snps, by = c("SNP" = "snpid")) %>%
  filter(Phenotype %in% mine)


# Construct GenomicRanges object ------------------------------------------
wranges <- with(wallace, GRanges(seqnames = chr,
                                 ranges = IRanges(start = pos, end = pos),
                                 pheno = trait,
                                 SNP = paste(chr, pos, sep = "_")))
sranges <- with(snps, GRanges(seqnames = chr,
                              ranges = IRanges(start = pos, end = pos),
                              Location = Location,
                              SNP = snpid))
stable <- as_data_frame(sranges@elementMetadata@listData)


# Control variables and parallel setup ------------------------------------
windows <- c(2e3, 5e3, 1e4, 1.5e4, 2e4, 3e4, 4e4, 5e4)
counts <- sig %>%
  filter(Type == t) %>%
  count(Phenotype, Location)
nperm <- 1000

cl <- makeCluster(8, outfile = "")
registerDoParallel(cl)


# Main loop ---------------------------------------------------------------
overlaps <- foreach(i = iter(windows), .packages = 'GenomicRanges') %dopar% {
  source("00.load-packages.R")
  cat("Window size:", i, "\n")
  
  wranges <- resize(wranges, width = 2*i + 1, fix = "center", use.names = FALSE)
  sranges <- resize(sranges, width = 2*i + 1, fix = "center", use.names = FALSE)
  w <- vector(mode = "numeric", length = nperm)
  
  for (j in 1:nperm) {
    cat(j, "\n")
    
    # Sample SNPs for each phenotype stratified by genomic location
    r <- vector(mode = "numeric")
    for (k in 1:nrow(counts)) {
      sample.snps <- stable %>%
        filter(Location == counts$Location[k]) %>%
        select(SNP) %>%
        unlist(use.names = FALSE)
      r <- c(r, sample(sample.snps, size = counts$n[k]))
    }
    
    ridx <- match(r, stable$SNP)
    rranges <- with(snps, GRanges(seqnames = chr[ridx],
                                  ranges = IRanges(start = pos[ridx], end = pos[ridx]),
                                  SNP = snpid[ridx],
                                  pheno = rep(counts$Phenotype, times = counts$n)))
    rranges <- resize(rranges, width = 2*i + 1, fix = "center", use.names = FALSE)
    
    # Check for overlaps
    rsnps <- findOverlaps(wranges, rranges, ignore.strand = TRUE)
    wind <- rsnps@from; rind <- rsnps@to
    d <- data_frame(MySNP = rranges@elementMetadata@listData$SNP[rind],
                    TheirSNP = wranges@elementMetadata@listData$SNP[wind],
                    MyPheno = rranges@elementMetadata@listData$pheno[rind],
                    TheirPheno = wranges@elementMetadata@listData$pheno[wind])
    w[j] <- with(d, sum(MyPheno == TheirPheno))
  }
  
  return(w)
}

names(overlaps) <- paste0("X", windows)
stopCluster(cl)

saveRDS(overlaps, paste0("data/overlap-permutations-iqr-", t, ".rds"))
