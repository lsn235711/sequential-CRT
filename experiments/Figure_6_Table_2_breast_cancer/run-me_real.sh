#!/bin/bash

for i in {1..100}
do
   sbatch --export=seed=$i script_real.sh
done


  
