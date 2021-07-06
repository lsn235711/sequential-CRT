#!/bin/bash

# Now run normal batch commands
seedXY=$1
seedX2=$2
ml R
Rscript sher_xvar.R $seedXY $seedX2
