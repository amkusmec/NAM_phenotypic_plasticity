### Natural log-transform the residuals for each phenotype.

source("00.load-packages.R")

# Data files and phenotype list -------------------------------------------
files <- list.files("./data/fwr-results", "*_IQR-fwr.rds")
phenos <- unlist(str_split(files, "-"))[seq(1, 2*length(files), 2)]
phenos <- unlist(str_split(phenos, "_"))[seq(1, 2*length(phenos), 2)]
traits <- readRDS("data/tidy_traitMatrix_IQR_AK.rds")

# Complete list of taxa because FarmCPU -----------------------------------
taxa <- readLines("data/genotypes_to_keep.txt")
taxa <- tibble(Taxa = taxa)

# Main loop ---------------------------------------------------------------
for (f in seq_along(files)) {
  cat(f, "\n")
  d <- readRDS(paste0("data/fwr-results/", files[f]))
  
  # Log-transform variances
  v <- tibble(Taxa = as.character(d$VAR), E = d$y - d$yhat[, 1]) %>%
    group_by(Taxa) %>%
    dplyr::summarise(VarE = log(var(E)))
  d2 <- tibble(Taxa = d$VARlevels, main = d$g[, 1], slope = d$b[, 1] + 1) %>%
    inner_join(., v, by = "Taxa") %>%
    right_join(., taxa, by = "Taxa")
  names(d2)[2:4] <- paste0(phenos[f], c("MainEffect", "Slope", "VarE"))
  
  saveRDS(d2, paste0("data/phenotypes/", phenos[f], "_IQR.rds"))
}
