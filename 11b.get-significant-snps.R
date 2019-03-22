setwd("~/gxe-gwas2/gwas-results")

library(tidyverse)
library(stringr)
library(qvalue)

threshold <- 0.01

sig <- list.files(".", pattern = "*.GWAS.Results.csv", full.names = TRUE,
                  recursive = TRUE) %>%
  map_df(function(x) {
    ty <- unlist(str_split(x, "/"))[2]
    ph <- unlist(str_split(x, "\\."))[4]
    read_csv(x, progress = FALSE) %>%
      mutate(Q.value = qvalue(p = P.value)$qvalues,
             Phenotype = gsub(ty, "", ph),
             Type = ty) %>%
      filter(Q.value <= threshold)
  })