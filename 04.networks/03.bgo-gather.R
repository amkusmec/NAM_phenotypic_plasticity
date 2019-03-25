### Gather the processed BINGO output and summarize it.

setwd("~/gxe-gwas2/data/bgo-edit")
source("../../00.load-packages.R")

# Biological process
files <- list.files("bp", "*.bgo", full.names = TRUE)
for (f in seq_along(files)) {
  stem <- unlist(str_split(gsub("bp/", "", files[f]), "-"))
  if (f == 1) {
    bp <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2])
  } else {
    bp <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2]) %>%
      bind_rows(bp, .)
  }
}

bp <- bp %>%
  select(Phenotype, Measure, GO.ID, Description, p.value:N) %>%
  mutate(Measure = gsub("MainEffect", "Mean Phenotype", Measure),
         Measure = gsub("Slope", "Linear Plasticity", Measure),
         Measure = gsub("VarE", "Non-linear Plasticity", Measure),
         Description = gsub(",", ";", Description))
write_csv(bp, "~/gxe-gwas2/supplementary/ppi_biological_process_enrichment.csv")

# Cellular component
files <- list.files("cc", "*.bgo", full.names = TRUE)
for (f in seq_along(files)) {
  stem <- unlist(str_split(gsub("cc/", "", files[f]), "-"))
  if (f == 1) {
    cc <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2])
  } else {
    cc <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2]) %>%
      bind_rows(cc, .)
  }
}

cc <- cc %>%
  select(Phenotype, Measure, GO.ID, Description, p.value:N) %>%
  mutate(Measure = gsub("MainEffect", "Mean Phenotype", Measure),
         Measure = gsub("Slope", "Linear Plasticity", Measure),
         Measure = gsub("VarE", "Non-linear Plasticity", Measure),
         Description = gsub(",", ";", Description))
write_csv(cc, "~/gxe-gwas2/supplementary/ppi_cellular_component_enrichment.csv")

# Molecular function
files <- list.files("mf", "*.bgo", full.names = TRUE)
for (f in seq_along(files)) {
  stem <- unlist(str_split(gsub("mf/", "", files[f]), "-"))
  if (f == 1) {
    mf <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2])
  } else {
    mf <- read_tsv(files[f]) %>%
      mutate(Phenotype = stem[1], Measure = stem[2]) %>%
      bind_rows(mf, .)
  }
}

mf <- mf %>%
  select(Phenotype, Measure, GO.ID, Description, p.value:N) %>%
  mutate(Measure = gsub("MainEffect", "Mean Phenotype", Measure),
         Measure = gsub("Slope", "Linear Plasticity", Measure),
         Measure = gsub("VarE", "Non-linear Plasticity", Measure),
         Description = gsub(",", ";", Description))
write_csv(mf, "~/gxe-gwas2/supplementary/ppi_molecular_function_enrichment.csv")
