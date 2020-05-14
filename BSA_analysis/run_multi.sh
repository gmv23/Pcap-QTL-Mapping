#!/usr/bin/bash

# Wrapper for running multipool on each chromosome and for different comparisons
# And then compiling output files together

path_to_multipool="/workdir/gmv23/bsa/multipool/multipool/mp_inference.py"

####### Run multipool in parallel ########

#Compare null: Rep 1 S vs Rep 2 S
mkdir output/s1_v_s2
seq -w 1 20 | parallel python $path_to_multipool -n 1056 data/chr{}_S1.txt data/chr{}_S2.txt -m contrast \
-c 100440 --plotFile output/s1_v_s2/chr{} -o output/s1_v_s2/chr{}.txt

#Compare null: Rep 1 T vs Rep 2 T
mkdir output/t1_v_t2
seq -w 1 20 | parallel $path_to_multipool -n 1056 data/chr{}_T1.txt data/chr{}_T2.txt -m contrast \
-c 100440 --plotFile output/t1_v_t2/chr{} -o output/t1_v_t2/chr{}.txt

#Compare contrast: Rep 1 S vs Rep 1 T
mkdir output/s1_v_t1
seq -w 1 20 | parallel python $path_to_multipool -n 1056 data/chr{}_S1.txt data/chr{}_T1.txt -m contrast \
-c 100440 --plotFile output/s1_v_t1/chr{} -o output/s1_v_t1/chr{}.txt

#Compare contrast: Rep 2 S vs Rep 2 T
mkdir output/s2_v_t2
seq -w 1 20 | parallel python $path_to_multipool -n 1056 data/chr{}_S2.txt data/chr{}_T2.txt -m contrast \
-c 100440 --plotFile output/s2_v_t2/chr{} -o output/s2_v_t2/chr{}.txt

########## Put it all together #########

if [ ! -d temp_joint ]; then
	mkdir temp_joint
fi

for i in $(seq -w 1 20)
	do
	#Get rid of header, get rid of second column, and sort each file
	tail -n +2 output/s1_v_s2/chr${i}.txt | cut -f1,3 | sort -k 1 -n > temp_joint/s1_v_s2_chr${i}.txt
	tail -n +2 output/t1_v_t2/chr${i}.txt | cut -f1,3 | sort -k 1 -n > temp_joint/t1_v_t2_chr${i}.txt
	tail -n +2 output/s1_v_t1/chr${i}.txt | cut -f1,3 | sort -k 1 -n > temp_joint/s1_v_t1_chr${i}.txt
	tail -n +2 output/s2_v_t2/chr${i}.txt | cut -f1,3 | sort -k 1 -n > temp_joint/s2_v_t2_chr${i}.txt

	echo "start joining"
	#Join files with different comparisons on the position field and add chromosome name
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,2.2' temp_joint/s1_v_s2_chr${i}.txt temp_joint/t1_v_t2_chr${i}.txt \
	| sort -k 1 -n | \
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,1.3,2.2' - temp_joint/s1_v_t1_chr${i}.txt \
	| sort -k 1 -n | \
	join -1 1 -2 1 -a 1 -a 2 -e 'NA' -o '0,1.2,1.3,1.4,2.2' - temp_joint/s2_v_t2_chr${i}.txt \
	| sort -k 1 -n | \
	awk -v i=$i '{print "Cp4.1LG"i"", $0}' > temp_joint/chr${i}_all.txt

done

#Now add header and concatenate all the files 
echo -e 'CHROM\tPOS\tS1_v_S2\tT1_v_T2\tS1_v_T1\tS2_v_T2' | cat - temp_joint/*_all.txt > multipool_lods.txt

#Remove temporary files
rm -r temp_joint
