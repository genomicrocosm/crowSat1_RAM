#!bin/bash/

### November 12th 2016 ###
### Matthias Weissensteiner ###
### mh.weissensteiner@gmail.com ###

### Created by Matthias Weissensteiner April 2016

# This script takes the output from the convert_to_input.sh script and combines
# it to a single file for the find_confs.sh step in LDhelmet

BASE="/proj/b2010059/nobackup/rho/SWE_extra"
inputfolder="${BASE}/raw_input"
outputfolder="${BASE}/whole_genome"

pops=SWE
cd $outputfolder
mkdir $pops
cd ${inputfolder}
for k in $(ls ${inputfolder})
do
echo $k
cd $k
	for i in $(ls ${pops}*INPUT)
		do
		echo $i
		sed -n '2~2p' $i > ${outputfolder}/${pops}/${i}.temp
		done
cd ..
done

cd $outputfolder/"$pops"
ls -1 *.temp | split -l 1000 -d - lists
for list in lists*; 
do paste -d "" $(cat $list) > merge${list##lists}; 
done
paste -d ""  merge* > seqs
sed -n '1~2p' ${inputfolder}/${k}/$i > IDs
paste IDs seqs | sed 's/\t/\n/g' > \
   ${outputfolder}/finished_files/${pops}.whole_genome.fasta
rm *temp
