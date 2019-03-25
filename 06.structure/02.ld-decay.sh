#!/bin/bash

### Calculate linkage disequilibrium between pairs of SNPs.

lib/plink --bfile data/geno_2M_all --r2 --ld-window-kb 300 --ld-window-r2 0 --out nam_ld
