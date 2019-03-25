### Run Bayesian Finlay-Wilkinson regression for all 23 phenotypes.
### This code assumes that you have 23 cores to run all regressions in parallel.

source("00.load-packages.R")
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Set-up control variables for parallel processing
ncores <- 23
cl <- makeCluster(ncores, outfile = "logs/bayesian_fwr_IQR.log")
registerDoParallel(cl)

# Load the trait matrix and split into a list
traitMatrix <- readRDS("data/tidy_traitMatrix_IQR_AK.rds")
phenos <- unique(traitMatrix$Phenotype)
traitList <- split(traitMatrix, traitMatrix$Phenotype)

### Bayesian FWR in parallel
# Control variables
burnin <- 1000
niter <- 51000

# Random seeds
seeds <- 47060859 + seq(0, ncores - 1)

fw <- foreach(ph = 1:23, .packages = 'FW') %dopar% {
    result <- with(traitList[[ph]], 
                   FW(y = Measure, VAR = Genotype, ENV = Environment, seed = seeds[ph],
                      method = "Gibbs", nIter = niter, burnIn = burnin, 
                      saveAt = paste0("data/gibbs-samples/", phenos[ph], "_IQR-Gibbs")))
    saveRDS(result, paste0("data/fwr-results/", phenos[ph], "_IQR-fwr.rds"))
}

stopCluster(cl)
