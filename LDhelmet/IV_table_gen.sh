#!/bin/bash

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### This script runs the table_gen step of LDhelmet to generate the likelihood
# lookup table with default parameters. 

pops=SWE

module load bioinfo-tools LDhelmet

BASE="/proj/b2010059/nobackup/rho/SWE_extra"
input_folder="${BASE}/conf_files"
output_folder="${BASE}/lk_tables"
theta=$(cat $input_folder/"$pops".theta)

time ldhelmet table_gen \
      --num_threads 8 \
      -t $theta \
      -r 0.0 0.1 10.0 1.0 100.0 \
      -c "$input_folder"/"$pops".conf \
      -o "$output_folder"/"$pops".lk
