setwd("~/gxe-gwas2")

source("00.load-packages.R")
library(corrplot)

gen_cor <- readRDS("gwas-results-iqr/genetic-correlations.rds")

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD",
                          "#4477AA"))

png("supplementary/gen_cor1.png", height = 12, width = 8, units = "in", res = 200)
par(mfrow = c(3, 2))
for (i in 1:6) {
  temp <- gen_cor[[i]]
  dimnames(temp) <- list(c("Mean", "Linear", "Non-linear"),
                         c("Mean", "Linear", "Non-linear"))
  corrplot(temp, method = "color", col = col(200), type = "upper",
           addCoef.col = "black", tl.col = "black", tl.cex = 1, tl.srt = 360,
           diag = FALSE, title = lbls[i], cl.pos = "b",
           mar = c(1.1, 1.1, 1.1, 1.1))
}
dev.off()

png("supplementary/gen_cor2.png", height = 12, width = 8, units = "in", res = 200)
par(mfrow = c(3, 2))
for (i in 7:12) {
  temp <- gen_cor[[i]]
  dimnames(temp) <- list(c("Mean", "Linear", "Non-linear"),
                         c("Mean", "Linear", "Non-linear"))
  corrplot(temp, method = "color", col = col(200), type = "upper",
           addCoef.col = "black", tl.col = "black", tl.cex = 1, tl.srt = 360,
           diag = FALSE, title = lbls[i], cl.pos = "b",
           mar = c(1.1, 1.1, 1.1, 1.1))
}
dev.off()

png("supplementary/gen_cor3.png", height = 12, width = 8, units = "in", res = 200)
par(mfrow = c(3, 2))
for (i in 13:18) {
  temp <- gen_cor[[i]]
  dimnames(temp) <- list(c("Mean", "Linear", "Non-linear"),
                         c("Mean", "Linear", "Non-linear"))
  corrplot(temp, method = "color", col = col(200), type = "upper",
           addCoef.col = "black", tl.col = "black", tl.cex = 1, tl.srt = 360,
           diag = FALSE, title = lbls[i], cl.pos = "b",
           mar = c(1.1, 1.1, 1.1, 1.1))
}
dev.off()

png("supplementary/gen_cor4.png", height = 12, width = 8, units = "in", res = 200)
par(mfrow = c(3, 2))
for (i in 19:23) {
  temp <- gen_cor[[i]]
  dimnames(temp) <- list(c("Mean", "Linear", "Non-linear"),
                         c("Mean", "Linear", "Non-linear"))
  corrplot(temp, method = "color", col = col(200), type = "upper",
           addCoef.col = "black", tl.col = "black", tl.cex = 1, tl.srt = 360,
           diag = FALSE, title = lbls[i], cl.pos = "b",
           mar = c(1.1, 1.1, 1.1, 1.1))
}
dev.off()
