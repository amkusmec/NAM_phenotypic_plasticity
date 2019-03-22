#!/bin/bash

lib/plink --bfile ~/gxe-gwas/data/geno_2M_all --r2 --ld-window-kb 300 --ld-window-r2 0 --out nam_ld
