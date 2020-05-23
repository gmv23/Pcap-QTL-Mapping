#!/usr/bin/bash

#Run all steps in bioinformatics analyses from whole-genome resequencing Squash BSA
#To turn raw reads into allele counts

set -e

SCRIPTS='/workdir/gmv23/squashQTL/Pcap-QTL-Mapping/BSA_bioinformatics'

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

if [ ! -d alignments ]; then
	mkdir alignments
	cd alignments
	cp "$SCRIPTS"/align.sh .
	mkdir bam_files
	bash align.sh ../trimming/trimmed_reads > align.log
	cd ..
fi

if [ ! -d variants ]; then
	mkdir variants
	cd variants
	mkdir bcf_files
	cp "$SCRIPTS"/call_variants.sh .
	bash call_variants.sh ../alignments/bam_files > call_variants.log
	cd ..
fi

if [ ! -d allele_counts ]; then
	mkdir allele_counts
	cd allele_counts
	cp "$SCRIPTS"/get_allele_counts.py .
	python get_allele_counts.py ../variants/final_geno_snps.recode.vcf allele_counts.txt
	cd ..
fi
