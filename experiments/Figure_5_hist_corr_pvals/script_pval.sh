#!/bin/bash

# Now run normal batch commands
seed=$1
ml R
Rscript sher_pval.R $seed
