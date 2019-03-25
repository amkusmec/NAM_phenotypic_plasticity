setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(car)

# Create tables for BC coefficients and normality tests -------------------
files <- list.files("./data/fwr-results", "*.rds")
phenos <- unlist(str_split(files, "-"))[seq(1, 2*length(files), 2)]

# Main loop ---------------------------------------------------------------
for (f in seq_along(files)) {
  cat(f, "\n")
  d <- readRDS(paste0("data/fwr-results/", files[f]))
  
  # QQ-plot of the residuals
  v <- tibble(Taxa = d$VAR, E = d$y - d$yhat[, 1])
  png(paste0("plots/qq-variance/residuals-", phenos[f], ".png"),
      width = 4, height = 4, units = "in", res = 150)
  qqPlot(v$E, envelope = 0.95, ylab = "Residuals", xlab = "Normal Quantiles",
         main = phenos[f])
  dev.off()
}
