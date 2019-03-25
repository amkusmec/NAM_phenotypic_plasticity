source("00.load-packages.R")

# Load and prep data ------------------------------------------------------
# Grab phenotypes from the trait matrix
traitMatrix <- readRDS("data/tidy_traitMatrix.rds")
phenos <- unique(traitMatrix$Phenotype)
rm(traitMatrix)

# Create tables to hold results
g.table <- tibble(Phenotype = "", Genotype = "", G = 0, Set = "")
b.table <- tibble(Phenotype = "", Genotype = "", B = 0, Set = "")
e.table <- tibble(Phenotype = "", Genotype = "", VE = 0, Set = "")


# Main loop ---------------------------------------------------------------
# Main loop to gather calculated values
for (i in seq_along(phenos)) {
  # Load the FW results
  new.res <- readRDS(paste0("data/fwr-results/", phenos[i], "-fwr.rds"))
  load(paste0("~/gxe-gwas/data/", phenos[i], "2-gibbs.RData"))
  
  # Calculate residual variances
  var.new <- tibble(Genotype = new.res$VAR, y = new.res$y, yhat = new.res$yhat[, 1]) %>%
    mutate(e = y - yhat) %>%
    group_by(Genotype) %>%
    dplyr::summarise(VE = var(e)) %>%
    mutate(Set = "New", Phenotype = phenos[i]) %>%
    select(Phenotype, Genotype, VE, Set)
  var.old <- tibble(Genotype = gibbs$VAR, y = gibbs$y, yhat = gibbs$yhat[, 1]) %>%
    mutate(e = y -yhat) %>%
    group_by(Genotype) %>%
    dplyr::summarise(VE = var(e)) %>%
    mutate(Set = "Old", Phenotype = phenos[i]) %>%
    select(Phenotype, Genotype, VE, Set)
  
  g.table <- bind_rows(list(g.table, 
                       tibble(Phenotype = phenos[i], Genotype = new.res$VARlevels, 
                              G = new.res$g[, 1], Set = "New"),
                       tibble(Phenotype = phenos[i], Genotype = gibbs$VARlevels, 
                              G = gibbs$g[, 1], Set = "Old")))
  b.table <- bind_rows(list(b.table,
                       tibble(Phenotype = phenos[i], Genotype = new.res$VARlevels, 
                              B = new.res$b[, 1] + 1, Set = "New"),
                       tibble(Phenotype = phenos[i], Genotype = gibbs$VARlevels,
                              B = gibbs$b[, 1] + 1, Set = "Old")))
  e.table <- bind_rows(list(e.table, var.new, var.old))
}

# Remove the initial value in each tibble
g.table <- g.table[-1, ]
b.table <- b.table[-1, ]
e.table <- e.table[-1, ]


# Calculate summary statistics --------------------------------------------
# Genotype main effect
inter1 <- g.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(Mean = mean(G)) %>%
  spread(key = Set, value = Mean) %>%
  dplyr::rename(New.Mean = New, Old.Mean = Old)
inter2 <- g.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(Var = var(G)) %>%
  spread(key = Set, value = Var) %>%
  dplyr::rename(New.Var = New, Old.Var = Old)
g.summary <- g.table %>%
  spread(key = Set, value = G) %>%
  group_by(Phenotype) %>%
  dplyr::summarise(Corr = cor(New, Old, method = "spearman", use = "pairwise")) %>%
  full_join(., inter1) %>%
  full_join(., inter2)

# Slope
inter1 <- b.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(Mean = mean(B)) %>%
  spread(key = Set, value = Mean) %>%
  dplyr::rename(New.Mean = New, Old.Mean = Old)
inter2 <- b.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(Var = var(B)) %>%
  spread(key = Set, value = Var) %>%
  dplyr::rename(New.Var = New, Old.Var = Old)
b.summary <- b.table %>%
  spread(key = Set, value = B) %>%
  group_by(Phenotype) %>%
  dplyr::summarise(Corr = cor(New, Old, method = "spearman", use = "pairwise")) %>%
  full_join(., inter1) %>%
  full_join(., inter2)

# Residual variance
# Will take 3 tables - median, 1st quartile, 3rd quartile
inter1 <- e.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(Median = median(VE)) %>%
  spread(key = Set, value = Median) %>%
  dplyr::rename(New.Median = New, Old.Median = Old)
inter2 <- e.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(X25 = quantile(VE, probs = 0.25)) %>%
  spread(key = Set, value = X25) %>%
  dplyr::rename(New.X25 = New, Old.X25 = Old)
inter3 <- e.table %>%
  group_by(Phenotype, Set) %>%
  dplyr::summarise(X75 = quantile(VE, probs = 0.75)) %>%
  spread(key = Set, value = X75) %>%
  dplyr::rename(New.X75 = New, Old.X75 = Old)
e.summary <- e.table %>%
  spread(key = Set, value = VE) %>%
  group_by(Phenotype) %>%
  dplyr::summarise(Corr = cor(New, Old, method = "spearman", use = "pairwise")) %>%
  full_join(., inter1) %>%
  full_join(., inter2) %>%
  full_join(., inter3)


# Visualize ---------------------------------------------------------------

b.table %>% group_by(Phenotype, Set) %>%
  dplyr::summarise(Mean = mean(B), SD = sd(B)) %>%
  ggplot(., aes(x = Set, y = Mean)) + theme_bw() + labs(y = "Slope") +
  geom_pointrange(aes(x = Set, ymin = Mean - SD, ymax = Mean + SD)) +
  facet_wrap("Phenotype", scales = "free_y")
ggsave("plots/slope-dispersion.png", width = 9, height = 9, units = "in", dpi = 200)

bind_rows(list(mutate(g.summary, Measure = "MainEffect"),
               mutate(b.summary, Measure = "Slope"),
               mutate(e.summary, Measure = "VarE"))) %>%
  ggplot(., aes(x = Corr, y = Phenotype)) + theme_bw() + xlim(-1, 1) +
  geom_point(size = I(3)) + labs(y = "", x = "Spearman Correlation") + facet_wrap("Measure")
ggsave("plots/new-old-correlations.png", width = 9, height = 8, units = "in", dpi = 200)

filter(b.table, Phenotype == "CobDiameter") %>%
  spread(key = Set, value = B) %>%
  ggplot(., aes(x = New, y = Old)) + theme_bw() + geom_hex() + 
  geom_abline(intercept = 0, slope = 1) +
  scale_fill_distiller(type = "seq", palette = "PuBuGn", direction = 1) +
  labs(x = "New Slope", y = "Old Slope", fill = "Count") + ggtitle("Cob Diameter")
ggsave("plots/CobDiameter-slope-comparison.png", width = 5, height = 5, units = "in", dpi = 200)

filter(b.table, Phenotype == "EarRowNumber") %>%
  spread(key = Set, value = B) %>%
  ggplot(., aes(x = New, y = Old)) + theme_bw() + geom_hex() + 
  geom_abline(intercept = 0, slope = 1) +
  scale_fill_distiller(type = "seq", palette = "PuBuGn", direction = 1) +
  labs(x = "New Slope", y = "Old Slope", fill = "Count") + ggtitle("Ear Row Number")
ggsave("plots/EarRowNumber-slope-comparison.png", width = 5, height = 5, units = "in", dpi = 200)

filter(b.table, Phenotype == "TasselPrimaryBranches") %>%
  spread(key = Set, value = B) %>%
  ggplot(., aes(x = New, y = Old)) + theme_bw() + geom_hex() + 
  geom_abline(intercept = 0, slope = 1) +
  scale_fill_distiller(type = "seq", palette = "PuBuGn", direction = 1) +
  labs(x = "New Slope", y = "Old Slope", fill = "Count") + ggtitle("Tassel Primary Branches")
ggsave("plots/TasselPrimaryBranches-slope-comparison.png", width = 5, height = 5, units = "in", dpi = 200)

filter(b.table, Phenotype == "UpperLeafAngle") %>%
  spread(key = Set, value = B) %>%
  ggplot(., aes(x = New, y = Old)) + theme_bw() + geom_hex() +
  geom_abline(intercept = 0, slope = 1) +
  scale_fill_distiller(type = "seq", palette = "PuBuGn", direction = 1) +
  labs(x = "New Slope", y = "Old Slope", fill = "Count") + ggtitle("Upper Leaf Angle")
ggsave("plots/UpperLeafAngle-slope-comparison.png", width = 5, height = 5, units = "in", dpi = 200)


# Count observations ------------------------------------------------------
traitMatrix <- readRDS("data/tidy_traitMatrix.rds")
oldMatrix <- readRDS("~/gxe-gwas/data/traitMatrix-nam.rds")
newCounts <- traitMatrix %>%
  group_by(Phenotype) %>%
  dplyr::summarise(N.New = n())
keyCounts <- oldMatrix %>%
  mutate(Key = paste0(Genotype, Phenotype)) %>%
  dplyr::count(Key) %>%
  filter(n >= 3)
oldCounts <- oldMatrix %>%
  mutate(Key = paste0(Genotype, Phenotype)) %>%
  filter(!is.na(Measure), 
         Genotype %in% traitMatrix$Genotype, 
         Key %in% keyCounts$Key) %>%
  select(-Key) %>%
  group_by(Phenotype) %>%
  dplyr::summarise(N.Old = n())
all.counts <- inner_join(newCounts, oldCounts)
