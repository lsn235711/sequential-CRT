#!/bin/bash

for i in {1..2000}
do
   sbatch --export=seed=$i script_hmm.sh
done


  
