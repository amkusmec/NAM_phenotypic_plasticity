setwd("~/gxe-gwas2/data")
library(bigmemory)
library(FarmCPUpp)


# Data input --------------------------------------------------------------
suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser()
parser$add_argument("-p", "--pheno", type = "character")
args <- parser$parse_args()
ph <- args$pheno

# Prep the phenotype data
cat(paste0("Loading phenotype data for ", args$pheno, "."), "\n")
myY <- readRDS(paste0("phenotypes/", args$pheno, "_IQR.rds"))

cat("Loading covariates.\n")
myQ <- readRDS("~/gxe-gwas/data/covariates6k.rds")
rownames(myQ) <- myQ$taxa
myQ <- myQ[, -1]

cat("Loading genotype map.\n")
myGM <- read.table("~/gxe-gwas/data/geno_2M.map", header = FALSE, 
                   stringsAsFactors = FALSE)
myGM <- myGM[, c(2, 1, 4)]
names(myGM) <- c("SNP", "Chromosome", "Position")
myGM$SNP <- paste0("X", myGM$SNP)

cat("Loading genotype matrix pointer.\n")
desc <- dget("geno_2M_pointer.desc")
myGD <- attach.big.matrix(desc)


# GWAS --------------------------------------------------------------------
for (j in 2:ncol(myY)) {
  setwd("~/gxe-gwas2/data")
  
  cat("FarmCPU on", names(myY)[j], "\n\n")
  start.time <- proc.time()[3]
  result <- farmcpu(Y = myY[, c(1, j)], GD = myGD, GM = myGM, CV = myQ, 
                    method.bin = "optimum", bin.size = c(5e3, 1e4, 5e4, 1e5), 
                    memo = "iqr", maxLoop = 20, ncores.glm = 32, ncores.reml = 16)
  end.time <- proc.time()[3]
  cat("Total runtime:", end.time - start.time, "(s).\n")
  
  if (j == 2) {
    setwd("~/gxe-gwas2/gwas-results-iqr/MainEffect")
  } else if (j == 3) {
    setwd("~/gxe-gwas2/gwas-results-iqr/Slope")
  } else {
    setwd("~/gxe-gwas2/gwas-results-iqr/VarE")
  }
  write_results(result)
  
  # Report and exit ---------------------------------------------------------
  runtime <- end.time - start.time
  system(paste0("python3 ~/mlmm-package/send_email.py -s '237 FarmCPU++ Finished - ", runtime, " s'"))
}
