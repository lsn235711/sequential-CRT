#!/bin/bash

# Now run normal batch commands
amp=$1
setting=$2
seed=$3
ml R
Rscript sher_comp.R $amp $setting $seed
