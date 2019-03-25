install.packages(c("tidyverse", "doParallel", "foreach", "iterators", 
                   "devtools", "car", "bigmemory", "biganalytics", 
                   "argparse", "Rcpp", "RcppEigen", "RcppParallel", 
                   "BiocManager", "VennDiagram", "scales", "lawstat", 
                   "readxl", "RColorBrewer", "scatterplot3d"))

devtools::install_github("lian0090/FW")
devtools::install_github("amkusmec/FarmCPUpp")
devtools::install_version("sommer", version = "2.9", 
                          repos = "http://cran.us.r-project.org")
BiocManager::install("qvalue", version = "3.8")
BiocManager::install("GenomicRanges", version = "3.8")
