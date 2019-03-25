setwd("~/gxe-gwas2")

source("00.load-packages.R")


# Prep candidate genes and orthologs --------------------------------------
candidates <- readRDS("gwas-results-iqr/candidate-genes.rds")
ortho <- read_csv("~/anno/at_orthologs.csv")
names(ortho) <- make.names(names(ortho))
ortho <- ortho %>% 
  filter(Arabidopsis.thaliana.gene.stable.ID != "") %>%
  inner_join(., candidates, by = c("Gene.stable.ID" = "GeneID"))


# Prep the GO annotations -------------------------------------------------
go <- read_delim("~/anno/ATH_GO_SLIM2.txt", delim = "\t", progress = FALSE) %>%
  select(Locus.name, GO.term, GO.ID, Ontology) %>%
  mutate(key = paste0(Locus.name, GO.ID)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-key)

go.translate <- go %>%
  group_by(GO.ID) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(GO.term, GO.ID)

go <- split(go, go$Ontology)
names(go) <- c("Cellular.Component", "Molecular.Function", "Biological.Process")


# Enrichment tests --------------------------------------------------------
enrichment <- vector(mode = "list")
for (o in names(go)) {
  cat("-", o, "\n")
  enrichment[[o]] <- vector(mode = "list")
  
  for (t in unique(ortho$Type)) {
    cat("--", t, "\n")
    
    # Get the selected genes
    myGenes <- ortho %>% filter(Type == t) %>%
      select(Arabidopsis.thaliana.gene.stable.ID) %>%
      unlist(use.names = FALSE) %>%
      unique(.)
    
    # Get the GO terms to test
    terms <- go[[o]] %>% filter(Locus.name %in% myGenes) %>%
      select(GO.ID) %>%
      unlist(use.names = FALSE) %>%
      unique(.)
    
    # Set up the structure to save results
    results <- matrix(0, nrow = length(terms), ncol = 6)
    
    for (j in seq_along(terms)) {
      term <- terms[j]
      
      white.draws <- nrow(filter(go[[o]], GO.ID == term, Locus.name %in% myGenes)) - 1
      white <- nrow(filter(go[[o]], GO.ID == term))
      black <- nrow(filter(go[[o]], GO.ID != term))
      draws <- length(myGenes)
      
      test <- phyper(white.draws, white, black, draws, lower.tail = FALSE)
      
      results[j, ] <- c(test, 0, white.draws, white, black, draws)
    }
    
    # Calculate corrected p-values
    results[, 2] <- p.adjust(results[, 1], method = "fdr")
    
    # Transform results into a data frame and add GO.term field
    colnames(results) <- c("p.value", "corr.p.value", "white.draw", "white", "black", "draws")
    results <- as_data_frame(results) %>%
      mutate(GO.ID = terms) %>%
      select(GO.ID, p.value:draws) %>%
      arrange(corr.p.value) %>%
      left_join(., go.translate, by = "GO.ID") %>%
      mutate(Type = t)
    
    enrichment[[o]][[t]] <- results
  }
}

# Save the results for further analysis
saveRDS(enrichment, "gwas-results-iqr/go-enrichment.rds")
