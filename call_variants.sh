#!/usr/bin/bash

# Sort and index bam files in parallel
# Then divide alignments by chromosome
# And run mpileup in parallel on each chromosome

#Provide path to directory containing bam files with no / at end
in_folder=$1

#set -x

# Sort input bam files in parallel
ls $in_folder/*.bam | parallel samtools sort {} -@ 4 -o {.}_sorted.bam

# Make file with chromosome names
grep ">" /workdir/gmv23/bsa/cpepo_genome/Cpepp_v4.1.chr.fa | sed 's/>//' > chromosome_names.txt

# Make file with names of input bam files
ls $in_folder/*sorted.bam > sorted_bams.txt

# Make file with names of bcf files for each chromosome
cat chromosome_names.txt | sed 's/\./_/g' \
sed 's/^/bcf_files\//' '| sed 's/$/.bcf/' > bcf_chrom_files.txt

# Index input bam files in parallel
parallel samtools index -@ 4 {} :::: sorted_bams.txt

# Run mpileup for each chromosome in parallel#
parallel --link bcftools mpileup -q 20 -b sorted_bams.txt -a AD,INFO/AD \
-f /workdir/gmv23/bsa/cpepo_genome/Cpepp_v4.1.chr.fa \
-r {1} -o {2} --threads 3 -O b :::: chromosome_names.txt \
:::: bcf_chrom_files.txt

# Now concatenate all the separate bcf files together
cd bcf_files
bcftools concat -f bcf_chrom_files.txt -O b --threads 24 -o final_allchrom.bcf

# Now call variants on pileup file 
bcftools call -c -v --threads 24 -O v -o final_geno.vcf final_allchrom.bcf

#Now filter to get two files:
#One with only SNPs for BSA
#One with Snps AND indels, filtered based on depth
