source("00.load-packages.R")

ld <- read_tsv("nam_ld.ld") %>%
  mutate(distance = (BP_B - BP_A)/1000,
         CHR_A = paste0("chr ", CHR_A),
         CHR_A = factor(CHR_A, levels = c("chr 1", "chr 2", "chr 3", "chr 4",
                                          "chr 5", "chr 6", "chr 7", "chr 8",
                                          "chr 9", "chr 10"), ordered = TRUE))

ggplot(ld, aes(x = distance, y = R2)) +
  geom_point(alpha = 0.3, shape = ".") +
  theme_bw() +
  facet_wrap("CHR_A", scales = "free_x", ncol = 2) +
  labs(x = "Distance (kb)", y = expression(r^2))
ggsave("supplementary/figs1-ld-decay.png", height = 8, width = 6.5, 
       units = "in", dpi = 300)
