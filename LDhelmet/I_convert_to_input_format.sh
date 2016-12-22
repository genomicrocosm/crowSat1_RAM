#!/bin/bash

### November 13th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

# Script to get fasta files per scaffold containing variant sites only from a 
# vcf file This script is part of the loop which is used to split up the 
# computational task and allocates one scaffold per job, the script used to 
# perform the loop is master-batchscript.sh


module load bioinfo-tools vcftools
module load plink

# the '$1' is to get the input from the loop, if many windows are to be analy
# zed, it is useful to split the bedfile up in smaller parts.
bedfile=$1 
BASE=/proj/b2010059/nobackup/rho/SWE_extra
inputfolder=/proj/b2010059/nobackup/Matthias_LD/input
outputfolder=${BASE}/raw_input
# phased vcf files, with one scaffold per vcf, are required as input
vcf_folder=${inputfolder}/phased_cleaned_vcf 
scaffold_appendix=list_genotypecrow.gatkHC.vqsrts99.phased.vcf
bedfolder=${inputfolder}/scaf_windows 

pop=SWE #population ID
cd $outputfolder
mkdir $bedfile

while read line;
 do
 # first split the bedfile up into scaffold and window coordinates, and
 # save them as variables
 scaffold=$(echo "$line" | cut -f 1 )
 from=$(echo "$line" | cut -f 2)
 to=$(echo "$line" | cut -f 3)
 echo $scaffold
 # then run vcftools with the --plink option, to the map and ped files from
 # the windows. the filtering options are minor allele frequency (0.1) and
 # missing data (no missing data allowed)
 vcf_file=${vcf_folder}/${scaffold}.CLEANED.${scaffold_appendix}
 vcftools --vcf ${vcf_file} \
 --chr ${scaffold} \
 --plink \
 --from-bp ${from} \
 --to-bp ${to} \
 --keep ${inputfolder}/pops_new/SWE \
 --out ${outputfolder}/${scaffold}/${pop}.${scaffold}.${from} \
 --maf 0.1 \ 
 --max-missing 1
 # then run plink to get the required 'fastphase' format
 plink --ped ${outputfolder}/${scaffold}/${pop}.${scaffold}.${from}.ped \
       --map ${outputfolder}/${scaffold}/${pop}.${scaffold}.${from}.map \
       --recode fastphase \
       --out ${outputfolder}/${scaffold}/${pop}.${scaffold}.${from}
 cd ${outputfolder}/${scaffold}
 # this requires then quite some bash acrobatics to get the required 
 # seq and snp pos file
 grep 'P' ${pop}.${scaffold}.${from}.chr-0.recode.phase.inp \
    | sed 's/ /\n/g' \
    | tail -n +2 > ${pop}.${scaffold}.${from}.snp.pos.file
 tail -n +4 ${pop}.${scaffold}.${from}.chr-0.recode.phase.inp \
    | sed 's/#/>/g' \
    | sed 's/ /_/g' > ${pop}.${scaffold}.${from}.temp 
 grep '>' ${pop}.${scaffold}.${from}.temp \
    > ${pop}.${scaffold}.${from}.ID_1
 awk '{print $1 "_2"}' ${pop}.${scaffold}.${from}.ID_1 \
    > ${pop}.${scaffold}.${from}.ID_2
 sed -n '2~3p' ${pop}.${scaffold}.${from}.temp \
    > ${pop}.${scaffold}.${from}.H1
 sed -n '3~3p' ${pop}.${scaffold}.${from}.temp \
    > ${pop}.${scaffold}.${from}.H2
 
 paste ${pop}.${scaffold}.${from}.ID_1 ${pop}.${scaffold}.${from}.H1 \
    | tr "\t" "\n" > ${pop}.${scaffold}.${from}.ID_1_H1
 paste ${pop}.${scaffold}.${from}.ID_2 ${pop}.${scaffold}.${from}.H2 \
    | tr "\t" "\n" > ${pop}.${scaffold}.${from}.ID_2_H2
 
 cat ${pop}.${scaffold}.${from}.ID_1_H1 \
     ${pop}.${scaffold}.${from}.ID_2_H2 \
     > ${pop}.${scaffold}.${from}.FINAL_INPUT
 rm *H1
 rm *H2
 rm *ID_1
 rm *ID_2
 rm *temp
 cd ..
done <${bedfolder}/${bedfile}


