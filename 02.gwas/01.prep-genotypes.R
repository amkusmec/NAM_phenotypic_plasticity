setwd("~/gxe-gwas2/data")

library(bigmemory)
library(biganalytics)

# Read the map first to get column names
myGM <- read.table("geno_2M.map", header = FALSE)
myGM <- myGM[, c(2, 1, 4)]
names(myGM) <- c("SNP", "Chromosome", "Position")
myGM$SNP <- paste0("X", myGM$SNP)

# Prep the genotype scores
myGD <- read.big.matrix("geno_2M.xmat", type = "double", sep = "\t", header = TRUE,
                        backingfile = "geno_2M.bin", descriptorfile = "geno_2M.desc",
                        extraCols = NULL, shared = TRUE, col.names = myGM$SNP,
                        ignore.row.names = FALSE, has.row.names = TRUE)

# Save the pointer for faster access when performing GWAS
desc <- describe(myGD)
dput(desc, "geno_2M_pointer.desc")
