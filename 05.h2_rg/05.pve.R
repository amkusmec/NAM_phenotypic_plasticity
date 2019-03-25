### Process LDAK output and summarise the estimated genomic heritabilities
### by annotation category.

source("00.load-packages.R")
library(scales)
library(RColorBrewer)


# Load the REML output ----------------------------------------------------
files <- list.files("data/reml", "*.reml", full.names = TRUE)
txt <- readLines(files[1])[11:16]
txt <- as.numeric(unlist(str_split(txt, " "))[seq(2, 3*length(txt), 3)])
stem <- unlist(str_split(files[1], "/"))[3]
d <- data_frame(Phenotype = rep(gsub("\\.reml", "", stem), 6),
                PVE = txt)

for (i in 2:length(files)) {
  txt <- readLines(files[i])[11:16]
  txt <- as.numeric(unlist(str_split(txt, " "))[seq(2, 3*length(txt), 3)])
  stem <- unlist(str_split(files[i], "/"))[3]
  d <- data_frame(Phenotype = rep(gsub("\\.reml", "", stem), 6),
                  PVE = txt) %>%
    bind_rows(d, .)
}


# Do some reformatting ----------------------------------------------------
categories <- c("Exon", "5kb Upstream", "5kb Downstream", "Intron",
                "Hypersensitivity Site", "Intergenic")
d$Category <- factor(rep(categories, times = 69), ordered = TRUE,
                     levels = categories)
d$Type <- ""
d[grepl("MainEffect", d$Phenotype), "Type"] <- "Mean Phenotype"
d[grepl("Slope", d$Phenotype), "Type"] <- "Linear Plasticity"
d[grepl("VarE", d$Phenotype), "Type"] <- "Non-linear Plasticity"
d$Type <- ordered(d$Type, 
                  levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"))
d <- d %>% mutate(Phenotype = rep(lbls, each = 18))


# Plots -------------------------------------------------------------------
ggplot(d, aes(x = Phenotype, y = PVE)) + theme_bw() +
  geom_bar(aes(fill = Category), stat = "identity") +
  labs(y = "Variance Explained", x = "", fill = "Annotation Category") +
  scale_fill_brewer(palette = "Set1") + 
  facet_wrap("Type", nrow = 3, strip.position = "left") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1))
ggsave("supplementary/figs5-pve.png", width = 8, height = 8, units = "in", dpi = 300)

d %>% group_by(Category, Type) %>%
  summarise(Mean = mean(PVE),
            StdErr = sd(PVE)/sqrt(length(lbls))) %>%
  ggplot(., aes(x = Category, y = Mean, colour = Type)) + theme_bw() +
  geom_pointrange(aes(ymin = Mean - StdErr, ymax = Mean + StdErr), 
                  position = position_dodge(width = 0.8)) +
  labs(y = "Mean Percent Variance Explained", x = "", fill = "") +
  scale_colour_manual(values = colors, breaks = levels(d$Type)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.3)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
ggsave("figures/fig2-mean-pve.pdf", width = 125, height = 142, units = "mm", 
       dpi = 300)
