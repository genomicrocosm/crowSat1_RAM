#!/bin/bash

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### This script is part of a wrapper script to split up the rjmcmc step
# of LDhelmet. As input it takes the likelihood table, the pade coefficients
# the mutation matrix, the sequence file and the snp position file. It is 
# run on default parameters with a block penalty of 50, a burn-in phase of
# 100000 and 1000000 iterations. 


module load bioinfo-tools LDhelmet

# population ID and scaffold are taken as input arguments from the 
# wrapper script.

pops=$1
scaffold=$2

BASE=/proj/b2010059/nobackup/rho/SWE_extra
lk_tables=${BASE}/lk_tables
pade_coefficients=${BASE}/pade_coefficients
mut_matrix=${BASE}/input
snp_seqs=${BASE}/raw_input/${scaffold}
snp_pos=${BASE}/raw_input/${scaffold}
output_folder=${BASE}/post_files

cd $output_folder
mkdir $scaffold
cd $snp_seqs
for i in $(ls "$pops".${scaffold}*.FINAL_INPUT)
	do
	input=$(echo $i | cut -f 1-3 -d ".")
	echo $input
	time ldhelmet rjmcmc \
	--num_threads 1  \
	-l "$lk_tables"/"$pops".lk  \
	-p "$pade_coefficients"/"$pops".pade  \
	-m "$mut_matrix"/zf_mut_mat.txt  \
	--snps_file "$snp_seqs"/"$input".FINAL_INPUT \
	--pos_file "$snp_pos"/"$input".snp.pos.file \
	-b 50 \
	--burn_in 100000  \
	-n 1000000  \
	-o "$output_folder"/"$scaffold"/"$input".post 
	done
