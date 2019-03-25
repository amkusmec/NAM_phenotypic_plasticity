setwd("~/gxe-gwas2")

source("00.load-packages.R")


# Prep the SNPs -----------------------------------------------------------
snps <- readRDS("data/snp-annotations.rds") %>%
  mutate(snpid = paste0("X", snpid))
names(snps)[4] = "Annotation"

sig <- readRDS("gwas-results-iqr/sig-snps.rds") %>%
  mutate(key = paste0(SNP, Type)) %>%
  group_by(key) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(-key) %>%
  left_join(., snps, by = c("SNP" = "snpid"))

msig <- sig %>% filter(Type == "MainEffect")
ssig <- sig %>% filter(Type == "Slope")
vsig <- sig %>% filter(Type == "VarE")


# Enrichment tests --------------------------------------------------------
# 1) Test significant SNPs against input distribution using 1D chi-square
probs <- snps %>% count(Annotation) %>% mutate(fract = n/sum(n))

mcat <- msig %>%
  count(Annotation) %>%
  mutate(fract = n/sum(n),
         totest = probs$fract*sum(n) >= 5,
         enrich = log2(fract/probs$fract),
         p.value = NA)
scat <- ssig %>%
  count(Annotation) %>%
  mutate(fract = n/sum(n),
         totest = probs$fract*sum(n) >= 5,
         enrich = log2(fract/probs$fract),
         p.value = NA)
vcat <- vsig %>%
  count(Annotation) %>%
  mutate(fract = n/sum(n),
         totest = probs$fract*sum(n) >= 5,
         enrich = log2(fract/probs$fract),
         p.value = NA)

with(mcat, chisq.test(n[totest], p = probs$n[totest], rescale.p = TRUE))
# 55.602, df = 4, p = 2.43e-11
with(scat, chisq.test(n[totest], p = probs$n[totest], rescale.p = TRUE))
# 24.81, df = 4, p = 5.493e-05
with(vcat, chisq.test(n[totest], p = probs$n[totest], rescale.p = TRUE))
# 7.6719, df = 3, p = 0.0533

# 2) Follow-up binomial tests
for (i in 1:nrow(probs)) {
  mcat$p.value[i] <- binom.test(mcat$n[i], n = sum(mcat$n), p = probs$fract[i], 
                                alternative = "two.sided")$p.value
  scat$p.value[i] <- binom.test(scat$n[i], n = sum(scat$n), p = probs$fract[i],
                                alternative = "two.sided")$p.value
  vcat$p.value[i] <- binom.test(vcat$n[i], n = sum(vcat$n), p = probs$fract[i],
                                alternative = "two.sided")$p.value
}

mcat <- mcat %>%
  filter(totest) %>%
  mutate(sig = p.value <= 0.05/n())
scat <- scat %>%
  filter(totest) %>%
  mutate(sig = p.value <= 0.05/n())
vcat <- vcat %>%
  filter(totest) %>%
  mutate(sig = p.value <= 0.05/n())

saveRDS(list(mcat, scat, vcat), "gwas-results-iqr/snp-enrichment.rds")
