#!/bin/bash -l

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

#### This is the script to run the find_confs step to build the haplotype 
# configuration file. As input it takes the whole genome fasta file from
# all the concatenated single window files. 



pops=SWE

module load bioinfo-tools LDhelmet
BASE="/proj/b2010059/nobackup/rho/SWE_extra"
input_folder="${BASE}/whole_genome/finished_files"
output_folder="${BASE}/conf_files"

time ldhelmet find_confs \
   --num_threads 8 \
   -w 50 \
   -o "$output_folder"/"$pops".conf "$input_folder"/"$pops".whole_genome.fasta 
