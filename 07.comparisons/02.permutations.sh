#!/bin/bash

Rscript 20.wallace-permutations.R --type MainEffect > logs/permute-MainEffect.log &
Rscript 20.wallace-permutations.R --type Slope > logs/permute-Slope.log &
Rscript 20.wallace-permutations.R --type VarE > logs/permute-VarE.log &
