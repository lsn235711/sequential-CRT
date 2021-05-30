#!/bin/bash

for i in {1..400}
do
   sbatch --export=seed=$i script_ar.sh
done


  
