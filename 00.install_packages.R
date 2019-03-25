install.packages(c("tidyverse", "doParallel", "foreach", "iterators", 
                   "devtools", "car", "bigmemory", "biganalytics", 
                   "argparse", "Rcpp", "RcppEigen", "RcppParallel", 
                   "BiocManager"))

devtools::install_github("lian0090/FW")
devtools::install_github("amkusmec/FarmCPUpp")
BiocManager::install("qvalue", version = "3.8")
BiocManager::install("GenomicRanges", version = "3.8")
