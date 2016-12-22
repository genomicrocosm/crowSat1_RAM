#!/bin/bash

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### This script is to run the pade step of LDhelmet to estimate the Padé 
# coefficients with default parameters. 


pops=SWE
module load bioinfo-tools LDhelmet

BASE="/proj/b2010059/nobackup/rho/SWE_extra"
input_folder="${BASE}/conf_files"
output_folder="${BASE}/pade_coefficients"
theta=$(cat $input_folder/"$pops".theta)

time ldhelmet pade \
       --num_threads 8 \
       -t $theta \
       -x 12 \
       --defect_threshold 40 \
       -c "$input_folder"/"$pops".conf \
       -o "$output_folder"/"$pops".pade
