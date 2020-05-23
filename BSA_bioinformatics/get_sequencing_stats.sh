#!usr/bin/bash

#This script is to get total #bp of
#raw reads
#trimmed and filtered reads
#reads aligning to genome
#Then the median SNP depth will accompany as well

path_raw_reads='/workdir/gmv23/bsa/data/round[12]'

path_filtered_reads=

for sample in D25 15_6015 1_S 1_T 2_S 2_T 2_R
do
	echo $sample
	gunzip -c $path_raw_reads/*Cpep_11_[A-Z]0[0-9]_$sample* | cat | \
	paste - - - - | cut -f 2 | tr -d '\n' | wc -c
done
