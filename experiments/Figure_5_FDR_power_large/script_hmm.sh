#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p candes,hns,normal,owners,pilanci
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=40960

# Now run normal batch commands
ml R
Rscript sher_hmm.R $seed
