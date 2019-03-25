source("00.load-packages.R")
library(scales)

sig <- readRDS("gwas-results-iqr/sig-snps.rds")

# Add total SNPs to the labels
totals <- sig %>% count(Phenotype)
lbls <- paste0(lbls, " (n=", totals$n, ")")

# Calculate phenotype-wise frequencies
percentages <- sig %>%
  group_by(Phenotype) %>%
  count(Type) %>%
  mutate(Total = sum(n), n = n/Total) %>%
  select(-Total) %>%
  ungroup()

# Reformat the type column
percentages <- percentages %>%
  mutate(Type = gsub("MainEffect", "Mean Phenotype", Type),
         Type = gsub("Slope", "Linear Plasticity", Type),
         Type = gsub("VarE", "Non-linear Plasticity", Type))
percentages$Type <- factor(percentages$Type, ordered = TRUE, 
                           levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"))

# Replace the phenotype column
totals <- percentages %>%
  count(Phenotype)
percentages <- percentages %>%
  mutate(Phenotype = rep(lbls, totals$nn))

# Plot
f3a <- ggplot(percentages, aes(x = Phenotype, y = n, fill = Type)) + theme_bw() +
  geom_bar(stat = "identity", position = position_stack(), colour = "black") +
  labs(x = "", y = "% Significant SNPs", fill = "") +
  geom_hline(yintercept = 0.5, linetype = 2, size = 2) +
  scale_fill_manual(values = colors) + #guides(fill = FALSE) +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 13, angle = 55, hjust = 1),
        legend.position = "top")
# ggsave("figures/measure-percentages.png", width = 10.5, height = 8, units = "in", dpi = 300)
