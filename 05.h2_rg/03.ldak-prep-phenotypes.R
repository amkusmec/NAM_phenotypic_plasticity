source("00.load-packages.R")


# Reformat phenotypes for LDAK --------------------------------------------
for (ph in phenos) {
  d <- readRDS(paste0("data/phenotypes/", ph, "_IQR.rds"))
  for (j in 2:ncol(d)) {
    write_tsv(d[, c(1, 1, j)], 
              path = paste0("data/phenotypes/", names(d)[j], ".txt"), 
              col_names = FALSE)
  }
}


# Reformat the principal components ---------------------------------------
readRDS("data/covariates6k.rds") %>%
  mutate(fid = taxa) %>%
  select(fid, taxa:PC3) %>%
  write_tsv(., path = "data/nam3pc.qcovar", col_names = FALSE)


# Write the analysis script -----------------------------------------------
txt <- c("#!/bin/bash", "")
counter <- 1

for (ph in phenos) {
  for (ty in c("MainEffect", "Slope", "VarE")) {
    if (counter %% 12 == 0) {
      txt <- c(txt,
               paste0("lib/ldak.4.9 --reml data/reml/", ph, ty, 
                      " --mgrm data/kinship-matrices/multi_grm.txt ",
                      "--kinship-details NO ",
                      "--covar data/nam3pc.qcovar ",
                      "--pheno data/phenotypes/", ph, ty, ".txt ",
                      "> logs/ldak-", ph, ty, ".log"))
    } else {
      txt <- c(txt,
               paste0("lib/ldak.4.9 --reml data/reml/", ph, ty, 
                      " --mgrm data/kinship-matrices/multi_grm.txt ",
                      "--kinship-details NO ",
                      "--covar data/nam3pc.qcovar ",
                      "--pheno data/phenotypes/", ph, ty, ".txt ",
                      "> logs/ldak-", ph, ty, ".log &"))
    }
    counter <- counter + 1
  }
}

writeLines(txt, "04.ldak.sh")
