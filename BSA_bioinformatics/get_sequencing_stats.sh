#!usr/bin/bash

#This script is to get total #bp of
#raw reads
#trimmed and filtered reads
#reads aligning to genome
#An R script will later append the median SNP read depth

#Note that this script is not parallelized at all and will take a very long time to run

path_raw_reads='/workdir/gmv23/bsa/data/round[12]'
path_filtered_reads='/workdir/gmv23/squashQTL/bsa_bioinformatics_2020-04-30/trimming/trimmed_reads'
path_bam_files='/workdir/gmv23/squashQTL/bsa_bioinformatics_2020-04-30/alignments/bam_files'

#Delete output file if it exists
if [ -f sequencing_stats.txt ]; then
	rm sequencing_stats.txt
fi

mkdir tmp

#Loop through different pools
for sample in D25 15_6015 1_S 1_T 2_S 2_T 2_R
do
	#Column for sample name
	echo -n $sample >> sequencing_stats.txt
	echo -n ',' >> sequencing_stats.txt

	#Column for raw read stats
	gunzip -c $path_raw_reads/*Cpep_11_[A-Z]0[0-9]_$sample* | cat | \
	paste - - - - | cut -f 2 | tr -d '\n' | wc -c | \
	xargs echo -n >> sequencing_stats.txt	
	echo -n ',' >> sequencing_stats.txt

	#Column for trimmed read stats
	gunzip -c $path_filtered_reads/*Cpep_11_[A-Z]0[0-9]_$sample* | cat | \
	paste - - - - | cut -f 2 | tr -d '\n' | wc -c | \
	xargs echo -n >> sequencing_stats.txt
	echo -n ',' >> sequencing_stats.txt

	#####Now get aligned bases
	#Add up bases in each subpool
	for file in $path_bam_files/[A-Z]0[0-9]_${sample}*sorted.bam
	do
		samtools view -b -q 20 $file | bedtools genomecov -ibam - -d | awk '{s+=$3}END{print s}' >> tmp/${sample}_aligned_bp.txt
	done
	#Add it all together and append to file
	awk '{s+=$0}END{print s}' tmp/${sample}_aligned_bp.txt | xargs echo >> sequencing_stats.txt

done

#rm -r tmp
