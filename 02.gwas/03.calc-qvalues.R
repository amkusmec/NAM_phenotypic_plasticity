setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(qvalue)

# for (ty in c("MainEffect", "Slope", "VarE")) {
#   setwd(paste0("~/gxe-gwas2/gwas-results/", ty))
#   files <- list.files(".", "*.GWAS.Results.csv")
#   for (f in files) {
#     cat(f, "\n")
#     res <- read_csv(f, progress = FALSE)
#     qobj <- qvalue(p = res$P.value)
#     res %>% mutate(Q.value = qobj$qvalues) %>%
#       write_csv(., path = f)
#   }
# }

list.files("gwas-results-iqr", "*GWAS.csv", full.names = TRUE, 
           recursive = TRUE) %>%
  map(function(f) {
    cat(f, "\n")
    res <- read_csv(f, progress = FALSE)
    qobj <- qvalue(p = res$p.value)
    res %>% mutate(q.value = qobj$qvalues) %>%
      write_csv(., path = f)
  })
