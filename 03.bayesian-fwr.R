setwd("~/gxe-gwas2")
# Master process ID: 51635

library(tibble)
library(doParallel)
library(foreach)
library(iterators)

# Set-up control variables for parallel processing
ncores <- 12
cl <- makeCluster(ncores, outfile = "logs/bayesian_fwr.log")
registerDoParallel(cl)

# Load the trait matrix and split into a list
traitMatrix <- readRDS("data/tidy_traitMatrix.rds")
phenos <- unique(traitMatrix$Phenotype)
traitList <- vector(mode = "list")
for (ph in phenos) {
  traitList[[ph]] <- subset(traitMatrix, Phenotype == ph)
}

### Bayesian FWR in parallel
# Control variables
burnin <- 5000
niter <- 50000

# Random seeds (chosen by running `round(runif(1, 1, 10^8))` in the console)
seeds <- 34884870 + seq(0, ncores - 1)

fw <- foreach(ph = iter(traitList), s = iter(seeds), .packages = 'FW') %dopar% {
  p <- unique(ph$Phenotype)
  with(ph, FW(y = Measure, VAR = Genotype, ENV = Environment, seed = s,
              method = "Gibbs", nIter = niter, burnIn = burnin, 
              saveAt = paste0("gibbs-samples/", p, "-Gibbs")))
}

stopCluster(cl)

# Save the results list
names(fw) <- phenos
saveRDS(fw, "data/fwr_results.rds")
