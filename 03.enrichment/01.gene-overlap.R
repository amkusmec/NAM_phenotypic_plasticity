source("00.load-packages.R")
library(VennDiagram)


# Load the candidate genes and remove duplicates --------------------------
d <- readRDS("gwas-results-iqr/candidate-genes.rds") %>%
  mutate(key = paste0(Type, GeneID)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-key)


# Pre-compute some quantities for the venn diagram ------------------------
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


# Venn diagram ------------------------------------------------------------
png("figures/gene-overlap.png", width = 5, height = 5, units = "in", res = 300)
draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12, n13 = n13, 
                 n23 = n23, n123 = n123, euler.d = FALSE, scaled = FALSE,
                 fill = colors, alpha = rep(0.5, 3), lwd = rep(3, 3), cex = rep(2, 7),
                 fontfamily = rep("sanserif", 7), cat.pos = c(345, 15, 180),
                 cat.dist = rep(0.05, 3), cat.cex = rep(1.5, 3),
                 cat.fontfamily = rep("sanserif", 3),
                 category = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"))
dev.off()


# Overlap enrichment tests ------------------------------------------------
# Load the full list of genes
gff <- read_delim("data/ZmB73_5b_FGS.gff", comment = "#", delim = "\t", 
                  progress = FALSE, col_names = FALSE) %>%
  filter(X3 == "gene", !is.na(X1)) %>%
  select(X1, X4, X5, X9) %>%
  mutate(X9 = unlist(str_split(X9, ";"))[seq(1, 3*n(), 3)]) %>%
  mutate(X9 = gsub('ID=', '', X9)) %>%
  arrange(X1, X4)
names(gff) <- c("chr", "start", "end", "gene")

# Test for enrichment of overlaps using Fisher's exact test
ms.test <- fisher.test(matrix(c(length(intersect(mgenes, sgenes)),
                                length(setdiff(sgenes, mgenes)),
                                length(setdiff(mgenes, sgenes)),
                                nrow(gff) - length(union(mgenes, sgenes))),
                              2, 2, byrow = TRUE))
mv.test <- fisher.test(matrix(c(length(intersect(mgenes, vgenes)),
                                length(setdiff(vgenes, mgenes)),
                                length(setdiff(mgenes, vgenes)),
                                nrow(gff) - length(union(mgenes, vgenes))),
                              2, 2, byrow = TRUE))
sv.test <- fisher.test(matrix(c(length(intersect(sgenes, vgenes)),
                                length(setdiff(vgenes, sgenes)),
                                length(setdiff(sgenes, vgenes)),
                                nrow(gff) - length(union(sgenes, vgenes))),
                              2, 2, byrow = TRUE))
tests <- data_frame(Group1 = c("Mean phenotype", "Mean phenotype", "Linear plasticity"),
                    Group2 = c("Linear plasticity", "Non-linear plasticity", "Non-linear plasticity"),
                    Odds.Ratio = c(ms.test$estimate, mv.test$estimate, sv.test$estimate),
                    Low95 = c(ms.test$conf.int[1], mv.test$conf.int[1], sv.test$conf.int[1]),
                    High95 = c(ms.test$conf.int[2], mv.test$conf.int[2], sv.test$conf.int[2]),
                    P.value = c(ms.test$p.value, mv.test$p.value, sv.test$p.value))
saveRDS(tests, "gwas-results-iqr/gene-overlap-enrichment.rds")
