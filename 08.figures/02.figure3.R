source("00.load-packages.R")
library(grid)
library(gridExtra)

source("03.enrichment/03.measure-by-phenotype.R")
source("07.comparisons/03.wallace-overlaps.R")

lay <- matrix(c(1, 1, 1, 2, 2), nrow = 1)
grob1 <- grobTree(ggplotGrob(f3a), 
                  textGrob("a", x = unit(0.03, "npc"), y = unit(0.975, "npc"), 
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
grob2 <- grobTree(ggplotGrob(f3b), 
                  textGrob("b", x = unit(0.03, "npc"), y = unit(0.975, "npc"),
                           hjust = "left", vjust = "top",
                           gp = gpar(fontface = "bold", fontsize = 14)))
g <- arrangeGrob(grob1, grob2, layout_matrix = lay)
ggsave("figures/fig3-measure-wallace.pdf", g, width = 300, height = 200, 
       units = "mm", dpi = 300)
