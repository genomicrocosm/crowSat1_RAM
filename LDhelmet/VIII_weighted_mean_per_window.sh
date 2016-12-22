#!/bin/bash

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### This script calculates the weighted mean of rho/bp per window with a size
# defined by a bedfile. 
# Input files are a file containing all merged rho estimates per window and
# a bedfile determining the window size. 


pops=$1
bed_file=$2
BASE=/proj/b2010059/nobackup/rho/SWE_extra
folder=${BASE}/concatenated_results

while read line; do
   scaffold=$(echo "$line" | cut -f 1 ) 
   from=$(echo "$line" | cut -f 2)
   to=$(echo "$line" | cut -f 3) 
   echo $from
   awk -v s=${scaffold} -v f=${from} -v t=${to} '$1 == s && $2 >= f' \
   ${folder}/${pops}.all_windows |\
   awk -v t=${to} '$3 <= t' > ${folder}/${pops}.${scaffold}.${from}.single_window
   start_window=$(tail -n +4 ${folder}/${pops}.${scaffold}.${from}.single_window  |\
     head -1 | awk '{print $2}') 
   end_window=$(tail -n +4 ${folder}/${pops}.${scaffold}.${from}.single_window |\
     tail -n 1 | awk '{print $3}')
   window_length=$(expr $end_window - $start_window)
   temp=$(awk '{print ($3-$2)*$4}' ${folder}/${pops}.${scaffold}.${from}.single_window | \
     awk '{ total += $1; count++ } END { print total }')	
   sum=$(printf "%.*f\n" 0 $temp)
   window=$(echo "$text_files"/"$scaffold"/$txt_file | cut -f 2 -d ".")
   echo $scaffold  $from $to $sum $window_length | \
     awk '{print $1 "\t" $2 "\t" $3 "\t"  $4/$5}' >> ${folder}/50kb/${pops}.${bed_file}_50kb_all

     #awk '{print $1 "\t" $2 "\t" $3 "\t"  $4/$5}' 
   rm ${folder}/${pops}.${scaffold}.${from}.single_window
done<${BASE}/scripts/pipeline/50kb/${bed_file}

