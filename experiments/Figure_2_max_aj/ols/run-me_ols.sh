#!/bin/bash

for i in {1..500}
do
   sbatch --export=seed=$i script_ols.sh
done


  
