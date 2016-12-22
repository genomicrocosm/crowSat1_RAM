#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# This is a wrapper which acts as a pipeline for getting the population rec-
# ombination rate rho from the program LDhelmet. It is adjusted for a 
# cluster which is operated with a SLURM system.
# written by Matthias Weissensteiner in January 2016
#------------------------------------------------------------------------------#

# Set up file paths:

#!/bin/bash

BASE=/proj/b2010059/nobackup/rho/SWE_extra
slurm_folder=${BASE}/scripts/pipeline/slurm_out
path_to_script=${BASE}/scripts/pipeline

#------------------------------------------------------------------------------#

# Stage 1 - build conf files

pops=SWE #population ID

sbatch -A b2014050 -p core -n 8 -t 40:00:00 \
    -J ${pops}.find_confs \
    -o $slurm_folder/find_confs_slurm/${pops}.find_confs.slurm \
    ${path_to_script}/III_find_confs.sh ${pops}

#------------------------------------------------------------------------------#

## Stage 2 - build likelihood tables

   sbatch -A b2014050 -p core -n 8 -t 150:00:00 \
    -J ${pops}.table_gen \
    -o ${slurm_folder}/lk_tables/${pops}.table_gen.slurm \
    $path_to_script/IV_table_gen.sh $pops

#------------------------------------------------------------------------------#

## Stage 3 - calculate Pad√ ®coefficients

   sbatch -A b2010059 -p core -n 16 -t 40:00:00 \
    -J "$pops".pade \
    -o $slurm_folder/$pops.pade.slurm \
    $path_to_script/V_pade.sh $pops

#------------------------------------------------------------------------------#

## Stage 4 - now calculate rho per window. 
# It is set up so that every scaffold is submitted that as a single job

for file in  ${BASE}/raw_input/*
    do
    scaffold=$(basename $file | cut -f 1 -d '.' )
    sbatch -A b2014050 -p core -n 1 -t 140:00:00 \
           -J ${pops}.${scaffold}.rjmcmc \
           -o ${slurm_folder}/rjmcmc/${pops}.${scaffold}.rjmcmc \
           "${path_to_script}/VI_rjmcmc.bash" "$pops" "$scaffold"
    done


#------------------------------------------------------------------------------#

## Stage 5 - This step converts the .post files into .text files

for folder in ${BASE}/post_files/* 
    do
    
    scaffold=$(basename "$folder")
    #echo $scaffold
    sbatch -A b2014050 -p core -n 1 -t 5:00:00 \
           -J ${scaffold}.post_to_text \
           -o ${slurm_folder}/post_to_text/${scaffold}.post_to_text \
           ${path_to_script}/VII_post_to_text.bash $scaffold
done


#------------------------------------------------------------------------------#
