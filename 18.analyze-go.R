setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(VennDiagram)

enrichment <- readRDS("gwas-results-iqr/go-enrichment.rds")
threshold <- 0.01

# List to store results
overlaps <- vector(mode = "list")

for (o in names(enrichment)) {
  # Get enriched terms at an FDR threshold
  m.enrich <- enrichment[[o]][["MainEffect"]] %>%
    filter(corr.p.value <= threshold)
  s.enrich <- enrichment[[o]][["Slope"]] %>%
    filter(corr.p.value <= threshold)
  v.enrich <- enrichment[[o]][["VarE"]] %>%
    filter(corr.p.value <= threshold)
  
  # 1) Generate vectors of overlapping and unique terms
  n1 <- setdiff(m.enrich$GO.term, union(s.enrich$GO.term, v.enrich$GO.term))
  n2 <- setdiff(s.enrich$GO.term, union(m.enrich$GO.term, v.enrich$GO.term))
  n3 <- setdiff(v.enrich$GO.term, union(m.enrich$GO.term, s.enrich$GO.term))
  n12 <- setdiff(intersect(m.enrich$GO.term, s.enrich$GO.term), v.enrich$GO.term)
  n13 <- setdiff(intersect(m.enrich$GO.term, v.enrich$GO.term), s.enrich$GO.term)
  n23 <- setdiff(intersect(s.enrich$GO.term, v.enrich$GO.term), m.enrich$GO.term)
  n123 <- intersect(m.enrich$GO.term, intersect(s.enrich$GO.term, v.enrich$GO.term))
  
  overlaps[[o]] <- list("M" = n1, "S" = n2, "V" = n3,
                        "MS" = n12, "MV" = n13, "SV" = n23,
                        "MSV" = n123)
  
  # 2) Make a Venn diagram of term overlaps
  n12 <- length(intersect(m.enrich$GO.ID, s.enrich$GO.ID))
  n13 <- length(intersect(m.enrich$GO.ID, v.enrich$GO.ID))
  n23 <- length(intersect(s.enrich$GO.ID, v.enrich$GO.ID))
  n123 <- length(intersect(m.enrich$GO.ID, intersect(s.enrich$GO.ID, v.enrich$GO.ID)))
  
  png(paste0("figures/", o, "-go-overlap.png"), width = 5, height = 5, units = "in", res = 300)
  draw.triple.venn(area1 = nrow(m.enrich), area2 = nrow(s.enrich), area3 = nrow(v.enrich),
                   n12 = n12, n13 = n13, n23 = n23, n123 = n123, 
                   euler.d = FALSE, scaled = FALSE, fill = colors, 
                   alpha = rep(0.5, 3), lwd = rep(3, 3), cex = rep(2, 7),
                   fontfamily = rep("sanserif", 7), cat.pos = c(345, 15, 180),
                   cat.dist = rep(0.05, 3), cat.cex = rep(1.5, 3),
                   cat.fontfamily = rep("sanserif", 3),
                   category = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"))
  dev.off()
  
  # 3) Construct a supplementary table
  enrich <- bind_rows(m.enrich, s.enrich, v.enrich) %>%
    select(Type, GO.ID, GO.term, p.value:draws) %>%
    mutate(Type = gsub("MainEffect", "Mean Phenotype", Type),
           Type = gsub("Slope", "Linear Plasticity", Type),
           Type = gsub("VarE", "Non-linear Plasticity", Type),
           GO.term = gsub(",", "\\.", GO.term))
  names(enrich)[1] <- "Measure"
  write_csv(enrich, path = paste0("supplementary/measure-", o, "-enrichment.csv"))
}

# Save the overlaps
saveRDS(overlaps, "gwas-results-iqr/go-term-overlaps.rds")
