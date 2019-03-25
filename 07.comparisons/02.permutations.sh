#!/bin/bash

### Run the permutation testing in parallel.

Rscript 07.comparisons/01.wallace-permutations.R --type MainEffect > logs/permute-MainEffect.log &
Rscript 07.comparisons/01.wallace-permutations.R --type Slope > logs/permute-Slope.log &
Rscript 07.comparisons/01.wallace-permutations.R --type VarE > logs/permute-VarE.log &
