### Assess overlaps between this project and Li et al.

library(tidyverse)
library(readxl)
library(GenomicRanges)

li <- read_excel("data/li_et_al_snps.xls", skip = 1)[-2502, ]
names(li) <- c("Population", "Trait", "SNP", "chr", "pos", "bpp")

window <- 1e5
nam <- filter(li, Population == "UA-NAM") %>%
  mutate(pos = pos*1e6)
da <- with(filter(nam, Trait == "DA"),
           GRanges(seqnames = chr,
                   ranges = IRanges(start = pos - window, end = pos + window)))
ds <- with(filter(nam, Trait == "DS"),
           GRanges(seqnames = chr,
                   ranges = IRanges(start = pos - window, end = pos + window)))
dacv <- with(filter(nam, Trait == "DACV"),
             GRanges(seqnames = chr,
                     ranges = IRanges(start = pos - window, end = pos + window)))
dscv <- with(filter(nam, Trait == "DSCV"),
             GRanges(seqnames = chr,
                     ranges = IRanges(start = pos - window, end = pos + window)))
asi <- with(filter(nam, Trait == "ASI"),
            GRanges(seqnames = chr,
                    ranges = IRanges(start = pos - window, end = pos + window)))

sig <- readRDS("gwas-results-iqr/sig-snps.rds") %>%
  select(Chromosome, Position, Phenotype, Type) %>%
  dplyr::rename(chr = Chromosome, pos = Position)
dts_me <- with(filter(sig, Phenotype == "DaystoSilk", Type == "MainEffect"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))
dts_sl <- with(filter(sig, Phenotype == "DaystoSilk", Type == "Slope"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))
dts_ve <- with(filter(sig, Phenotype == "DaystoSilk", Type == "VarE"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))
dtt_me <- with(filter(sig, Phenotype == "DaysToTassel", Type == "MainEffect"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))
dtt_sl <- with(filter(sig, Phenotype == "DaysToTassel", Type == "Slope"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))
dtt_ve <- with(filter(sig, Phenotype == "DaysToTassel", Type == "VarE"),
               GRanges(seqnames = chr,
                       ranges = IRanges(start = pos - window, end = pos + window)))

r_da <- reduce(da)
r_dacv <- reduce(dacv)
r_ds <- reduce(ds)
r_dscv <- reduce(dscv)

length(findOverlaps(dacv, da))
length(findOverlaps(dscv, ds))

length(findOverlaps(dtt_me, da))
length(findOverlaps(dtt_sl, da))
length(findOverlaps(dtt_ve, da))

length(findOverlaps(dtt_me, dacv))
length(findOverlaps(dtt_sl, dacv))
length(findOverlaps(dtt_ve, dacv))

length(findOverlaps(dts_me, ds))
length(findOverlaps(dts_sl, ds))
length(findOverlaps(dts_ve, ds))

length(findOverlaps(dts_me, dscv))
length(findOverlaps(dts_sl, dscv))
length(findOverlaps(dts_ve, dscv))

# GEMMA results for DTS
library(qvalue)
# No significant hits for Main Effect?
gsl <- read_tsv("gemma/output/DaystoSilk_IQR_Slope.assoc.txt") %>%
  mutate(q = qvalue(p_lrt)$qvalues) %>%
  filter(q <= 0.01)
gsl_range <- with(gsl,
                  GRanges(seqnames = chr,
                          ranges = IRanges(start = ps - window, end = ps + window)))
# No significant hits for Residual Variance
gsl_range <- reduce(gsl_range)
length(findOverlaps(gsl_range, ds))
length(findOverlaps(gsl_range, dscv))
length(findOverlaps(gsl_range, dts_me))
length(findOverlaps(gsl_range, dts_sl))
length(findOverlaps(gsl_range, dts_ve))
