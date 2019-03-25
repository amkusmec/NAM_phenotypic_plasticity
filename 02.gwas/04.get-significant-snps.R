### Idenfity significant SNPs at FDR <= 1%.

setwd("~/gxe-gwas2/gwas-results-iqr")

source("../00.load-packages.R")

threshold <- 0.01
sig <- list.files(".", pattern = "*.GWAS.csv", full.names = TRUE,
                  recursive = TRUE) %>%
  map_df(function(x) {
    ty <- unlist(str_split(x, "/"))[2]
    ph <- unlist(str_split(x, "\\."))[3]
    read_csv(x, progress = FALSE) %>%
      mutate(Phenotype = gsub(ty, "", ph),
             Type = ty) %>%
      filter(q.value <= threshold)
  })
saveRDS(sig, "sig-snps.rds")

# Save a version as a supplementary table
sig <- sig %>%
  dplyr::select(Phenotype, Type, SNP:Position, estimate, stderr, p.value, q.value) %>%
  mutate(Type = gsub("MainEffect", "Mean Phenotype", Type),
         Type = gsub("Slope", "Linear Plasticity", Type),
         Type = gsub("VarE", "Non-linear Plasticity", Type))
names(sig)[2] <- "Measure"
write_csv(sig, "../supplementary/significant-snps.csv")
