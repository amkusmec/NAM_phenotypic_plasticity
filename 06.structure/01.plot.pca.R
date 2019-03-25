### Plot the first 3 principal components of the genotype matrix. The NAM
### family structure is clearly captured by these components.

library(stringr)
library(scatterplot3d)

pca <- readRDS("data/covariates6k.rds")
pca <- subset(pca, !grepl("M", pca$taxa))
pca$fid <- unlist(str_split(pca$taxa, "E"))[seq(1, 2*nrow(pca), 2)]
pca$fid <- factor(pca$fid)

x <- colors(distinct = TRUE)[-1]
idx <- seq(1, length(x), floor(length(x)/26))
x <- x[idx[-27]]

png("supplementary/figs3-pca-scatterplot.png", height = 6, width = 6, units = "in", 
    res = 300)
par(mar = c(4.1, 4.1, 1.1, 1.1), mfrow = c(2, 2))
with(pca, plot(PC1, PC2, pch = 19, cex = 0.5, col = x[fid],
               xlab = "PC1 (2.6%)", ylab = "PC2 (2.1%)"))
with(pca, plot(PC1, PC3, pch = 19, cex = 0.5, col = x[fid],
               xlab = "PC1 (2.6%)", ylab = "PC3 (2.1%)"))
with(pca, plot(PC2, PC3, pch = 19, cex = 0.5, col = x[fid],
               xlab = "PC2 (2.1%)", ylab = "PC3 (2.1%)"))
with(pca, scatterplot3d(PC2, PC3, PC1, color = x[fid], mar = c(3.1, 2.1, 1.1, 2.1), 
                        pch = 19, cex.symbols = 0.5, cex.axis = 0.5,
                        las = 2, xlab = "PC2 (2.1%)", ylab = "PC3 (2.1%)",
                        zlab = "PC1 (2.6%)"))
dev.off()

