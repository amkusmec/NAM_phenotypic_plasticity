### Assemble main text Figure 4.

source("00.load-packages.R")
library(grid)
library(gridExtra)
library(VennDiagram)

source("03.enrichment/02.window-overlaps.R")


# Gene overlaps -----------------------------------------------------------
d <- readRDS("gwas-results-iqr/candidate-genes.rds") %>%
  mutate(key = paste0(Type, GeneID)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-key)

mgenes <- d %>% filter(Type == "MainEffect") %>% select(GeneID) %>% unlist(use.names = FALSE)
sgenes <- d %>% filter(Type == "Slope") %>% select(GeneID) %>% unlist(use.names = FALSE)
vgenes <- d %>% filter(Type == "VarE") %>% select(GeneID) %>% unlist(use.names = FALSE)

n1 <- length(mgenes)
n2 <- length(sgenes)
n3 <- length(vgenes)
n12 <- length(intersect(mgenes, sgenes))
n13 <- length(intersect(mgenes, vgenes))
n23 <- length(intersect(sgenes, vgenes))
n123 <- length(intersect(mgenes, intersect(sgenes, vgenes)))

gene.venn <- draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12,
                              n13 = n13, n23 = n23, n123 = n123, euler.d = FALSE,
                              scaled = FALSE, fill = colors,
                              lwd = rep(3, 3), cex = rep(2, 7),
                              fontfamily = rep("Helvetica", 7), 
                              cat.pos = c(345, 15, 180),
                              cat.dist = rep(0.05, 3), cat.cex = rep(1.5, 3),
                              cat.fontfamily = rep("Helvetica", 3),
                              category = c("Mean Phenotype", "Linear Plasticity",
                                           "Non-linear plasticity"), ind = FALSE)


# GO overlaps -------------------------------------------------------------
enrichment <- readRDS("gwas-results-iqr/go-enrichment.rds")
m.enrich <- enrichment$Biological.Process$MainEffect %>%
  filter(corr.p.value <= 0.01)
s.enrich <- enrichment$Biological.Process$Slope %>%
  filter(corr.p.value <= 0.01)
v.enrich <- enrichment$Biological.Process$VarE %>%
  filter(corr.p.value <= 0.01)

n1 <- length(m.enrich$GO.ID)
n2 <- length(s.enrich$GO.ID)
n3 <- length(v.enrich$GO.ID)
n12 <- length(intersect(m.enrich$GO.ID, s.enrich$GO.ID))
n13 <- length(intersect(m.enrich$GO.ID, v.enrich$GO.ID))
n23 <- length(intersect(s.enrich$GO.ID, v.enrich$GO.ID))
n123 <- length(intersect(m.enrich$GO.ID, intersect(s.enrich$GO.ID, v.enrich$GO.ID)))

go.venn <- draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12,
                            n13 = n13, n23 = n23, n123 = n123, euler.d = FALSE,
                            scaled = FALSE, fill = colors,
                            lwd = rep(3, 3), cex = rep(2, 7),
                            fontfamily = rep("Helvetica", 7), 
                            cat.pos = c(345, 15, 180),
                            cat.dist = rep(0.05, 3), cat.cex = rep(1.5, 3),
                            cat.fontfamily = rep("Helvetica", 3),
                            category = c("Mean Phenotype", "Linear Plasticity",
                                         "Non-linear plasticity"), ind = FALSE)


# Assemble the figure -----------------------------------------------------
lay <- matrix(c(1, 2, 3), nrow = 1)
grob1 <- grobTree(gene.venn,
                  textGrob("a", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)),
                  textGrob("*", x = unit(0.55, "npc"), y = unit(0.785, "npc"),
                           gp = gpar(fontsize = 16)))
grob2 <- grobTree(ggplotGrob(f4b),
                  textGrob("b", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
grob3 <- grobTree(go.venn,
                  textGrob("c", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
g <- arrangeGrob(grob1, grob2, grob3, layout_matrix = lay)
ggsave("figures/fig4-overlaps.pdf", g, width = 400, height = 150,
       units = "mm", dpi = 300)
