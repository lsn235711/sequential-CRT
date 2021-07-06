#!/bin/bash

# Slurm parameters
PART=candes,hns,normal,owners,pilanci  # Partition names
MEMO=40960                     # Memory required (10G)
TIME=24:00:00               # Time required (1h)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -p "$PART" --time="$TIME

LOGS=logs
mkdir -p $LOGS

for seed in {1..500}; do
    # Script to be run
    SCRIPT="script_pval.sh $seed"
    # Define job name
    JOBN="pval_"$seed
    OUTF=$LOGS"/"$JOBN".out"
    ERRF=$LOGS"/"$JOBN".err"
    # Assemble slurm order for this job
    ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
    # Print order
    echo $ORD
    # Submit order
    $ORD
    #./$SCRIPT
done



  
