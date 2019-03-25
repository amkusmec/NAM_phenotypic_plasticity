source("00.load-packages.R")


# Function for constructing regular expressions ---------------------------
get.pattern <- function(query) {
  if (grepl("GRMZM", query)) {
    pat <- paste0(query, "_P")
  } else {
    query <- unlist(str_split(query, "G"))[1]
    pat <- paste0(query, "GP")
  }
  
  return(pat)
}


# Prep data ---------------------------------------------------------------
candidates <- readRDS("gwas-results/candidate-genes.rds")
net <- read_delim("~/anno/Zm_ppi.sif", col_names = FALSE, delim = "\t")
names(net) <- c("Node1", "Type", "Node2")


# Construct subnetworks ---------------------------------------------------
phenos <- unique(candidates$Phenotype)
types <- unique(candidates$Type)

for (i in seq_along(phenos)) {
  cat(phenos[i], "\n")
  
  for (k in seq_along(types)) {
    genes <- candidates %>%
      filter(Phenotype == phenos[i], Type == types[k]) %>%
      select(GeneID) %>% unlist(use.names = FALSE)
    
    if (length(genes) > 0) {
      pat <- get.pattern(genes[1])
      temp <- net[grepl(pat, net$Node1) | grepl(pat, net$Node2), ]
      
      for (j in 2:length(genes)) {
        pat <- get.pattern(genes[j])
        temp <- bind_rows(temp, net[grepl(pat, net$Node1) | grepl(pat, net$Node2), ])
      }
      
      write_tsv(temp, path = paste0("data/networks/", phenos[i], "-", types[k], ".sif"),
                col_names = FALSE)
    }
  }
}


# Construct an aggregate network ------------------------------------------
files <- list.files("data/networks", "*.sif", full.names = TRUE)
agg <- read_delim(files[1], delim = "\t", col_names = FALSE) %>%
  mutate(Phenotype = gsub("\\.sif", "", unlist(str_split(files[1], "/"))[3]))
for (f in 2:length(files)) {
  agg <- read_delim(files[f], delim = "\t", col_names = FALSE) %>%
    mutate(Phenotype = gsub("\\.sif", "", unlist(str_split(files[f], "/"))[3])) %>%
    bind_rows(., agg)
}


# Extract lists of genes for each phenotype -------------------------------
phenos <- unique(agg$Phenotype)
batches <- vector(mode = "character")
for (p in phenos) {
  genes <- c(subset(agg, Phenotype == p)$X1, subset(agg, Phenotype == p)$X3)
  genes <- unique(genes)
  batches <- c(batches, p, genes, "batch")
}


# Filter the aggregate network --------------------------------------------
agg <- agg %>%
  mutate(key = paste0(X1, X2)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(X1:X3)


# Output ------------------------------------------------------------------
write_tsv(agg, path = "data/networks/aggregate.sif", col_names = FALSE)
write_lines(batches, path = "data/networks/batches.txt")
