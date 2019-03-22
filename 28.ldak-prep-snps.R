setwd("~/gxe-gwas2")

source("00.load-packages.R")

anno <- readRDS("data/snp-annotations.rds")
for (ty in levels(anno$type)) {
  anno %>% filter(type == ty) %>%
    select(snpid) %>% unlist(use.names = FALSE) %>%
    writeLines(., paste0("data/kinship-matrices/", gsub(" ", "", ty), "-snps.txt"))
}
