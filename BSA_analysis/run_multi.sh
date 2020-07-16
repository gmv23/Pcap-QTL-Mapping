#!/usr/bin/bash

# Wrapper for running multipool on each chromosome and for different comparisons
# And then compiling output files together

path_to_multipool="/workdir/gmv23/bsa/multipool/multipool/mp_inference.py"

###################################### Run multipool in parallel ######################

#Write wrapper function to save stand out from multipool to file
#Arguments:
# 1) path to multipool
# 2) input file 1
# 3) input file 2
# 4) output base

launch_multipool(){
python2 $1 -n 1056 $2 $3 -m contrast -c 125636 --plotFile $4 -o $4.txt &> $4.log
}
export -f launch_multipool

#Compare null: Rep 1 S vs Rep 2 S
mkdir output/s1_v_s2
seq -w 1 20 | parallel launch_multipool $path_to_multipool data/s1_v_s2/chr{}_S1.txt \
data/s1_v_s2/chr{}_S2.txt output/s1_v_s2/chr{}

#Compare null: Rep 1 T vs Rep 2 T
mkdir output/t1_v_t2
seq -w 1 20 | parallel launch_multipool $path_to_multipool data/t1_v_t2/chr{}_T1.txt \
data/t1_v_t2/chr{}_T2.txt output/t1_v_t2/chr{}

#Compare contrast: Rep 1 S vs Rep 1 T
mkdir output/s1_v_t1
seq -w 1 20 | parallel launch_multipool $path_to_multipool data/s1_v_t1/chr{}_S1.txt \
data/s1_v_t1/chr{}_T1.txt output/s1_v_t1/chr{}

#Compare contrast: Rep 2 S vs Rep 2 T
mkdir output/s2_v_t2
seq -w 1 20 | parallel launch_multipool $path_to_multipool data/s2_v_t2/chr{}_S2.txt \
data/s2_v_t2/chr{}_T2.txt output/s2_v_t2/chr{}

########## Put together results with LOD scores across genome from output txt files #########

if [ ! -d temp_joint ]; then
	mkdir temp_joint
fi

for i in $(seq -w 1 20)
	do
	#Get rid of header, get rid of second column, and sort each file
	tail -n +2 output/s1_v_s2/chr${i}.txt | cut -f1,3 | sort -k 1 > temp_joint/s1_v_s2_chr${i}.txt
	tail -n +2 output/t1_v_t2/chr${i}.txt | cut -f1,3 | sort -k 1 > temp_joint/t1_v_t2_chr${i}.txt
	tail -n +2 output/s1_v_t1/chr${i}.txt | cut -f1,3 | sort -k 1 > temp_joint/s1_v_t1_chr${i}.txt
	tail -n +2 output/s2_v_t2/chr${i}.txt | cut -f1,3 | sort -k 1 > temp_joint/s2_v_t2_chr${i}.txt

	echo "start joining"
	#Join files with different comparisons on the position field and add chromosome name
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,2.2' temp_joint/s1_v_s2_chr${i}.txt temp_joint/t1_v_t2_chr${i}.txt \
	| sort -k 1 | \
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,1.3,2.2' - temp_joint/s1_v_t1_chr${i}.txt \
	| sort -k 1 | \
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,1.3,1.4,2.2' - temp_joint/s2_v_t2_chr${i}.txt \
	| sort -k 1 -n | \
	awk -v i=$i '{print "Cp4.1LG"i"", $0}' > temp_joint/chr${i}_all.txt

done

#Now add header and concatenate all the files 
echo -e 'CHROM\tPOS\tS1_v_S2\tT1_v_T2\tS1_v_T1\tS2_v_T2' | cat - temp_joint/*_all.txt > multipool_lods.txt

#Remove temporary files
rm -r temp_joint

######################## Pull out peak locations and CIs from log files ###############

#Write header line
echo -e 'CHROM,Rep1_start,Rep1_end,Rep1_peak,Rep2_start,Rep2_end,Rep2_peak' >> CIs_and_peaks.txt

#Loop through chromosomes
for chrom in $(seq -w 1 20)
do
	#First print chomosome
	echo -n 'Cp4.1LG' >> CIs_and_peaks.txt
	echo -n $chrom >> CIs_and_peaks.txt
	echo -ne ',' >> CIs_and_peaks.txt

	#Pull CI start and stop and peak location for rep 1
	grep '90% credible interval' output/s1_v_t1/chr$chrom.log | \
	sed -r 's/^.*?spans ([0-9]*?) ([0-9]*?) .*?mode: ([0-9]*)$/\1,\2,\3,/' | \
	xargs echo -n >> CIs_and_peaks.txt

	#Pull CI start and stop and peak location for rep 2
	grep '90% credible interval' output/s2_v_t2/chr$chrom.log | \
	sed -r 's/^.*?spans ([0-9]*?) ([0-9]*?).*?mode: ([0-9]*)$/\1,\2,\3/' | \
	xargs echo >> CIs_and_peaks.txt
done
