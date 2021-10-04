# The Sequential CRT
This repository contains the code to reproduce the numerical results in our paper: **Deploying the Conditional Randomization Test
in High Multiplicity Problems**.

## Folders
* sequential_CRT/: contains the main functions that implement the sequential CRT proposed in our paper. 
* examples/: contains examples that demonstrate the usage of the functions.
* experiments/: contains the code to reproduce the numerical results in our paper.

## Replicating the experiments
* For Figure 1 and 6, the ".R" files in the corresponding folders can be directly executed on a laptop.
* For the other figures and tables, one needs to run the "run-me\*.sh" file on a computing cluster, and then run the "make_plot.R" or "make_table.R" file to make the plot or table. The "script\*.sh" file is specifically designed for the [Stanford Sherlock](https://www.sherlock.stanford.edu/) cluster. One may need to modify the files correspondingly. 
