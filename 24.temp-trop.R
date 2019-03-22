setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(lawstat)


# Control variables and storage structures --------------------------------
phenos <- c("DaystoSilk", "DaysToTassel", "GDDDaystoSilk", "GDDDaystoTassel")
traits <- c("Days to Silk", "Days to Tassel", "GDD to Silk", "GDD to Tassel")
dispersion <- vector(mode = "list")
subgroups <- data_frame(fid = c("Z001", "Z002", "Z003", "Z004", "Z005", "Z006", "Z007", "Z008",
                                "Z009", "Z010", "Z011", "Z012", "Z013", "Z014", "Z015", "Z016",
                                "Z017", "Z018", "Z019", "Z020", "Z021", "Z022", "Z023", "Z024",
                                "Z025", "Z026"),
                        tt = c("Temperate", "Tropical", "Tropical", "Temperate", "Tropical",
                               "Tropical", "Tropical", "Tropical", "Tropical", "Temperate",
                               "Temperate", "Tropical", "Tropical", "Temperate", "Mixed", "Mixed",
                               "Temperate", "Tropical", "Temperate", "Tropical", "Tropical",
                               "Temperate", "Temperate", "Temperate", "Mixed", "Tropical"))
subgroups$tt <- factor(subgroups$tt)


# Function for calculating dispersion -------------------------------------
qcd <- function(x) {
  q <- quantile(x)
  (q[4] - q[2])/(q[4] + q[2])
}


# Get dispersion ----------------------------------------------------------
for (i in seq_along(phenos)) {
  gibbs <- readRDS(paste0("data/fwr-results/", phenos[i], "-fwr.rds"))
  
  # Calculate dispersion family-wise
  dispersion[[phenos[i]]] <- data_frame(taxa = gibbs$VARlevels, Slope = gibbs$b[, 1] + 1) %>%
    mutate(fid = unlist(str_split(taxa, "E"))[seq(1, 2*n(), 2)]) %>%
    group_by(fid) %>%
    summarise(Disp = qcd(Slope)) %>%
    left_join(., subgroups, by = "fid") %>%
    mutate(trait = traits[i])
}


# Tests -------------------------------------------------------------------
# 1) Levene's test for homogeneity of variances
# 2) Kruskal-Wallis test for among-group differences
for (i in seq_along(phenos)) {
  cat(traits[i], "\n")
  print(levene.test(dispersion[[phenos[i]]]$Disp, dispersion[[phenos[i]]]$tt, 
                    location = "median"))
  print(kruskal.test(Disp ~ tt, data = dispersion[[i]]))
}


# Follow-up tests ---------------------------------------------------------
# GDDDaystoSilk
silk <- data_frame(Group1 = c("Mixed", "Mixed", "Temperate"),
                   Group2 = c("Temperate", "Tropical", "Tropical"),
                   W = NA, p.value = NA)
for (i in 1:nrow(silk)) {
  z <- wilcox.test(subset(dispersion[[3]], tt == silk$Group1[i])$Disp,
                   subset(dispersion[[3]], tt == silk$Group2[i])$Disp)
  silk$W[i] <- z$statistic
  silk$p.value[i] <- z$p.value
}

# GDDDaystoTassel
tass <- data_frame(Group1 = c("Mixed", "Mixed", "Temperate"),
                   Group2 = c("Temperate", "Tropical", "Tropical"),
                   W = NA, p.value = NA)
for (i in 1:nrow(tass)) {
  z <- wilcox.test(subset(dispersion[[4]], tt == tass$Group1[i])$Disp,
                   subset(dispersion[[4]], tt == tass$Group2[i])$Disp)
  tass$W[i] <- z$statistic
  tass$p.value[i] <- z$p.value
}


# Plot --------------------------------------------------------------------
segments <- tibble(x = c(1, 2, 1, 2),
                   y = c(0.27, 0.3, 0.15, 0.2),
                   xend = c(2, 3, 2, 3),
                   yend = c(0.27, 0.30, 0.15, 0.2),
                   trait = c("GDD to Silk", "GDD to Silk", "GDD to Tassel", "GDD to Tassel"))
stars <- tibble(t = c("*", "*", "*", "*"), x = c(1.5, 2.5, 1.5, 2.5),
                y = c(0.29, 0.32, 0.17, 0.22),
                trait = c("GDD to Silk", "GDD to Silk", "GDD to Tassel", "GDD to Tassel"))
bind_rows(dispersion[[1]], dispersion[[2]], dispersion[[3]], dispersion[[4]]) %>%
  ggplot(., aes(tt, Disp)) + ylim(0, 0.35) +
  geom_point(size = 2) + theme_bw() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data = segments, size = 1) +
  geom_text(aes(x = x, y = y, label = t), data = stars, size = 7) +
  labs(y = "Quartile Coefficient of Determination", x = "") +
  facet_wrap("trait")
ggsave("supplementary/figs4-flowering-disp.png", width = 8, height = 8, units = "in", dpi = 600)
