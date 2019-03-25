source("00.load-packages.R")
library(gridExtra)

grobs1 <- 1:6 %>% map(function(i) {
  test <- gen_cor[[i]]
  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(test[, 1], test[, 2], test[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(lbls[i])
})

png("supplementary/figs2a-gen-cor.png", width = 8, height = 12, units = "in", res = 300)
grid.arrange(grobs = grobs1, ncol = 2)
dev.off()

grobs2 <- 7:12 %>% map(function(i) {
  test <- gen_cor[[i]]
  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(test[, 1], test[, 2], test[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(lbls[i])
})

png("supplementary/figs2b-gen-cor.png", width = 8, height = 12, units = "in", res = 300)
grid.arrange(grobs = grobs2, ncol = 2)
dev.off()

grobs3 <- 13:18 %>% map(function(i) {
  test <- gen_cor[[i]]
  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(test[, 1], test[, 2], test[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(lbls[i])
})

png("supplementary/figs2c-gen-cor.png", width = 8, height = 12, units = "in", res = 300)
grid.arrange(grobs = grobs3, ncol = 2)
dev.off()

grobs4 <- 19:23 %>% map(function(i) {
  test <- gen_cor[[i]]
  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(test[, 1], test[, 2], test[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(lbls[i])
})

png("supplementary/figs2d-gen-cor.png", width = 8, height = 12, units = "in", res = 300)
grid.arrange(grobs = grobs4, ncol = 2)
dev.off()
