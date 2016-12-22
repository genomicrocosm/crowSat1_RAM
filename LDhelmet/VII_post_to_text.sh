#!/bin/bash

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### This script converts the post files from the rjmcmc step into txt files.

module load bioinfo-tools LDhelmet

scaffold=$1

BASE=/proj/b2010059/nobackup/rho/SWE_extra
post_files=${BASE}/post_files
text_files=${BASE}/text_files

cd ${text_files}
mkdir $scaffold
cd ${post_files}/${scaffold}
for i in $(ls); do
    ldhelmet post_to_text -m -p 0.025 -p 0.50 -p 0.975 \
    -o ${text_files}/${scaffold}/${i}.txt \
    ${post_files}/${scaffold}/$i
    done

