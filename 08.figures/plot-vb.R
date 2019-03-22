setwd("~/gxe-gwas2")

library(tidyverse)
library(stringr)

files <- list.files("data/fwr-results", "*-fwr.rds", full.names = TRUE)[seq(2, 46, 2)]
vb <- files %>%
  map_df(function(f) {
    gibbs <- read_rds(f)
    data_frame(NumEnv = length(gibbs$ENVlevels),
               Lower = gibbs$var_b - gibbs$SD.var_b,
               Var = gibbs$var_b,
               Upper = gibbs$var_b + gibbs$SD.var_b)
  }) %>%
  do.call("rbind", .) %>%
  t() %>%
  as_data_frame() %>%
  mutate(Phenotype = unlist(str_split(files, "/"))[seq(3, 3*length(files), 3)],
         Phenotype = unlist(str_split(Phenotype, "-"))[seq(1, 2*length(files), 2)])

ggplot(vb, aes(x = NumEnv, y = Var)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) +
  labs(x = "Number of Environments", y = expression(sigma[b]^2)) +
  theme_bw()
ggsave("supplementary/env-vs-sigmab.png", width = 5, height = 5, units = "in", 
       dpi = 300)

ggplot(vb, aes(x = Phenotype, y = Var)) +
  geom_pointrange(aes(ymin = Lower, ymax = Upper, colour = factor(NumEnv))) +
  scale_colour_brewer(type = "seq", palette = 7) +
  labs(x = "Phenotype", y = expression(sigma[b]^2), colour = "# Environments") +
  coord_flip() +
  theme_bw()
ggsave("supplementary/sigmab-by-phenotype.png", width = 5, height = 8, 
       units = "in", dpi = 300)

cor.test(vb$NumEnv, vb$Var)
cor.test(vb$NumEnv[-18], vb$Var[-18])
