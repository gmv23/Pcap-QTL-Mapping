#!/usr/bin/bash

#Run all steps in bioinformatics analyses from whole-genome resequencing Squash BSA
#To turn raw reads into allele counts
#Script modified from run_all_steps.sh in sexual_pop scripts

set -e

SCRIPTS='/workdir/gmv23/squashQTL/Pcap-QTL-Mapping'

if [ -z $1 ]; then
        DATE=$(date +%F)
else
    	DATE=$1
        echo "Date received"
fi

if [ ! -d bsa_bioinformatics_$DATE ]; then
        mkdir bsa_bioinformatics_$DATE
fi

cd bsa_bioinformatics_$DATE

#Copy current version of this script to directory
cp "$SCRIPTS"/reads_to_allele_counts.sh .

if [ ! -d trimming ]; then
        mkdir trimming
        cd trimming

	#Make file with paths to all raw reads from both rounds 1 and 2 of sequencing
	ls /workdir/gmv23/bsa/data/round1/*fastq.gz >> read_paths.txt
	ls /workdir/gmv23/bsa/data/round2/*fastq.gz >> read_paths.txt

	#Run fastp in parallel on all reads
	cp "$SCRIPTS"/trim.sh .
	mkdir trimmed_reads

	bash trim.sh read_paths.txt trimmed_reads > trimming.log
	cd ..
fi

