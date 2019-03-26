## Distinct genetic architectures for phenotype means and plasticities in *Zea mays*

### Citation

Kusmec, A., S. Srinivasan, D. Nettleton, and P.S. Schnable. (2017) Distinct genetic architectures for phenotype means and plasticities in *Zea mays*. *Nature Plants*, **3**: 715-23. (https://www.nature.com/articles/s41477-017-0007-7)

You may also be interested in [this](https://www.nature.com/articles/s41477-017-0012-x) accompanying perspective.

### Introduction

This project uses data on 23 phenotypes measured on ~5,000 maize inbred lines in 4-11 environments (location-years) to study phenotypic plasticity. We found that the candidate genes for average phenotype across environments and for response to environments are different in *Zea mays*.

### Description

Prior to running any of the scripts in the repository, please run `00.install_packages.R` to install all R packages required to reproduce these analyses. The necessary directory structure has been preserved through the use of `.gitignore` files. For the raw data files, please see the paper for links.

1. **phenotypes**: These scripts process raw phenotypic data and perform Bayesian Finlay-Wilkinson regression to quantify phenotypic plasticity.
2. **gwas**: These scripts perform GWAS for phenotypic plasticity, apply multiple testing correction, and generate lists of significant genetic markers and candidate genes.
3. **enrichment**: These scripts summarize the GWAS results, phenotype distributions, and perform enrichment for gene ontology annotations.
4. **networks**: These scripts examine gene ontology enrichment results based on a *Zea mays* protein-protein interaction network. These scripts will need to be supplemented with analyses in Cytoscape as detailed in the manuscript.
5. **h2_rg**: These scripts calculate genomic heritability based on annotations of genetic markers and compute additive genetic heritabilities using a multi-trait mixed model.
6. **structure**: These scripts summarize population structure and the rate of linkage disequilibrium decay in the population under study.
7. **comparisons**: These scripts conduct comparisons between the results of this study and two other studies, one of which explicitly studied phenotypic variability across environments.
8. **figures**: These scripts create composite main text and supplementary figures.

### License

This repository is free and open source for use in research and is licensed under the terms of [GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/#).
