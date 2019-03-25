source("00.load-packages.R")
library(qvalue)

list.files("gwas-results-iqr", "*GWAS\\.csv", full.names = TRUE, 
           recursive = TRUE) %>%
  map(function(f) {
    cat(f, "\n")
    res <- read_csv(f, progress = FALSE)
    qobj <- qvalue(p = res$p.value)
    res %>% mutate(q.value = qobj$qvalues) %>%
      write_csv(., path = f)
  })
