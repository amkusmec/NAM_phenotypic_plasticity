source("00.load-packages.R")
library(FarmCPUpp)
library(qvalue)

files.m <- list.files("gwas-results-iqr/MainEffect", "*GWAS.csv", full.names = TRUE)
files.s <- list.files("gwas-results-iqr/Slope", "*GWAS.csv", full.names = TRUE)
files.v <- list.files("gwas-results-iqr/VarE", "*GWAS.csv", full.names = TRUE)

for (i in 1:length(files.m)) {
  cat(i, "\n")
  
  res.m <- read_csv(files.m[i])
  names(res.m)[4] <- "p.value"
  qobj <- qvalue(res.m$p.value)
  thresh.m <- max(qobj$pvalues[qobj$qvalues <= 0.01], na.rm = TRUE)
  
  res.s <- read_csv(files.s[i])
  names(res.s)[4] <- "p.value"
  qobj <- qvalue(res.s$p.value)
  thresh.s <- max(qobj$pvalues[qobj$qvalues <= 0.01], na.rm = TRUE)
  
  res.v <- read_csv(files.v[i])
  names(res.v)[4] <- "p.value"
  qobj <- qvalue(res.v$p.value)
  thresh.v <- max(qobj$pvalues[qobj$qvalues <= 0.01], na.rm = TRUE)
  
  png(paste0("supplementary/figs6", letters[i], "-", phenos[i], ".png"), width = 9,
      height = 5.5, units = "in", res = 300)
  layout(matrix(c(1, 1, 2, 3, 3, 4, 5, 5, 6), nrow = 3, byrow = TRUE))
  
  manhattan_plot(res.m, colors = c("black", "grey60"), cutoff = thresh.m,
                 main = "Mean Phenotype", pch = 19, cex = 0.5)
  qq_plot(res.m, pch = 19, cex = 0.6)
  
  manhattan_plot(res.s, colors = c("black", "grey60"), cutoff = thresh.s,
                 main = "Linear Plasticity", pch = 19, cex = 0.5)
  qq_plot(res.s, pch = 19, cex = 0.6)
  
  manhattan_plot(res.v, colors = c("black", "grey60"), cutoff = thresh.v,
                 main = "Non-linear Plasticity", pch = 19, cex = 0.5)
  qq_plot(res.v, pch = 19, cex = 0.6)
  
  mtext(paste0("  ", letters[i], ". ", lbls[i]), side = 3, line = -2, 
        outer = TRUE, font = 2, adj = 0)
  
  dev.off()
}
