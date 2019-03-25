### Calculate quartile coefficient of dispersion for the plasticity measures.

source("00.load-packages.R")


# Set-up storage variables ------------------------------------------------
nenv <- vector(mode = "numeric", length = length(lbls))
slope <- data_frame(Phenotype = lbls, Disp = 0, Median = 0, Mean = 0, Var = 0)
vare <- data_frame(Phenotype = lbls, Disp = 0, Median = 0, Mean = 0, Var = 0)


# Function to calculate the quantile coefficient of dispersion ------------
qcd <- function(x) {
  q <- quantile(x, na.rm = TRUE)
  (q[4] - q[2])/(q[4] + q[2])
}


# Calculate QCD for each phenotype ----------------------------------------
for (i in seq_along(phenos)) {
  gibbs <- readRDS(paste0("data/fwr-results/", phenos[i], "_IQR-fwr.rds"))
  nenv[i] <- length(gibbs$ENVlevels)
  
  slope$Disp[i] <- qcd(gibbs$b[, 1] + 1)
  slope$Median[i] <- median(gibbs$b[, 1] + 1)
  slope$Mean[i] <- mean(gibbs$b[, 1] + 1)
  slope$Var[i] <- var(gibbs$b[, 1] + 1)
  
  v <- data_frame(taxa = gibbs$VAR, e = gibbs$y - gibbs$yhat[, 1]) %>%
    group_by(taxa) %>%
    summarise(VarE = var(e))
  vare$Disp[i] <- qcd(v$VarE)
  vare$Median[i] <- median(v$VarE)
  vare$Mean[i] <- mean(v$VarE)
  vare$Var[i] <- var(v$VarE)
}


# Plot --------------------------------------------------------------------
disp <- bind_rows(slope, vare) %>%
  mutate(Type = factor(rep(c("Linear Plasticity", "Non-linear Plasticity"), 
                           each = length(lbls)), ordered = TRUE),
         Phenotype = paste0(Phenotype, " (n=", rep(nenv, times = 2), ")"),
         Phenotype = factor(Phenotype, ordered = TRUE, 
                            levels = unique(Phenotype)))

ggplot(disp, aes(x = Phenotype, y = Disp, colour = Type, shape = Type)) +
  geom_point(size = I(4)) + theme_bw() + ylim(0, 1) +
  scale_colour_manual(values = colors[-1]) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Quartile Coefficient of Dispersion", x = "", colour = "", shape = "")
ggsave("figures/fig1-dispersion.pdf", width = 180, height = 200, 
       units = "mm", dpi = 300)
