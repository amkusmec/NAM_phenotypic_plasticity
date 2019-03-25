install.packages(c("tidyverse", "doParallel", "foreach", "iterators", 
                   "devtools", "car", "bigmemory", "biganalytics", 
                   "argparse", "Rcpp", "RcppEigen", "RcppParallel", 
                   "BiocManager", "VennDiagram", "scales", "lawstat", 
                   "readxl", "RColorBrewer", "scatterplot3d"))

### Software for Bayesian Finlay-Wilkinson regresssion
devtools::install_github("lian0090/FW")

### Parallelized, C++ version of FarmCPU for efficiency
devtools::install_github("amkusmec/FarmCPUpp")

### SOMMER has undergone major revisions since the publication of this paper.
### Use of the latest version will require major revision of the code to 
### calculate additive genetic correlations.
devtools::install_version("sommer", version = "2.9", 
                          repos = "http://cran.us.r-project.org")

BiocManager::install("qvalue", version = "3.8")
BiocManager::install("GenomicRanges", version = "3.8")
