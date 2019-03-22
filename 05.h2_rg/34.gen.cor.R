setwd("~/gxe-gwas2")

library(tidyverse)
library(sommer)

# Prepare the kinship matrix
A.mat <- read.table("~/gxe-gwas/data/geno_9k_kin.txt", header = FALSE, skip = 1)
taxa <- A.mat[, 1]
A.mat <- as.matrix(A.mat[, -1])
rownames(A.mat) <- taxa; colnames(A.mat) <- taxa

gen_cor <- list.files("data/phenotypes", "*_IQR.rds", full.names = TRUE) %>%
  map(function(f) {
    Y <- readRDS(f)
    Y <- Y[!is.na(Y[, 2]), ]
    A <- A.mat[rownames(A.mat) %in% Y$Taxa, colnames(A.mat) %in% Y$Taxa]
    Y <- Y[, -1]
    
    Za <- diag(nrow(Y))
    ETA.A <- list(add = list(Z = Za, K = A))
    ans.A <- mmer(Y = Y, Z = ETA.A, MVM = TRUE, method = "EMMA")
    
    gvc <- ans.A$var.comp$Vu
    sd.gvc <- as.matrix(sqrt(diag(gvc)))
    prod.sd <- sd.gvc %*% t(sd.gvc)
    gvc/prod.sd
  })

saveRDS(gen_cor, "gwas-results-iqr/genetic-correlations.rds")



Y <- readRDS("data/phenotypes/DaystoSilk_IQR.rds")
Y <- Y[!is.na(Y[, 2]),]
me <- Y$DaystoSilkMainEffect
sl <- Y$DaystoSilkSlope
ve <- Y$DaystoSilkVarE
taxa <- Y$Taxa
idx <- match(taxa, rownames(A.mat), nomatch = 0)
A.mat <- A.mat[idx, idx]
ans.A <- mmer2(cbind(DaystoSilkMainEffect, DaystoSilkSlope, DaystoSilkVarE) ~ 1, 
               random = ~ g(Taxa), G = list(Taxa = A.mat), data = Y, MVM = TRUE)
