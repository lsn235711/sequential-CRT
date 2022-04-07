#!/bin/bash

# Now run normal batch commands
B=$1
setting=$2
seed=$3
amp=$4
ml R
Rscript sher_comp.R $B $setting $seed $amp
