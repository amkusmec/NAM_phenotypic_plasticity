#!/bin/bash

# Construct annotation-specific kinship matrices using GCTA

cd ~/gxe-gwas2/data/kinship-matrices

../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract 5kbdownstream-snps.txt --make-grm-bin --out 5kbdownstream > ../../logs/kinship-5kbdownstream.log &
../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract 5kbupstream-snps.txt --make-grm-bin --out 5kbupstream > ../../logs/kinship-5kbupstream.log &
../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract Exon-snps.txt --make-grm-bin --out exon > ../../logs/kinship-exon.log &
../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract Hypersensitivitysite-snps.txt --make-grm-bin --out hypersensitivity > ../../logs/kinship-hs.log &
../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract Intergenic-snps.txt --make-grm-bin --out intergenic > ../../logs/kinship-intergenic.log &
../../lib/gcta/gcta64 --bfile ../geno_2M_all --extract Intron-snps.txt --make-grm-bin --out intron > ../../logs/kinship-intron.log &
