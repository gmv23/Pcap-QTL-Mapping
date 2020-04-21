#!/usr/bin/bash

#Align reads to genome in parallel

in_folder=$1
out_folder=$2

# Make new files with base names and file paths for input to bwa in parallel
if [ -f $out_folder/base_names.txt ]; then
	rm $out_folder/base_names.txt
fi
if [ -f $out_folder/file_names.txt ]; then
	rm $out_folder/file_names.txt
fi

# Get base names from round1 files (flow cell HWYTYBGX5) (Kind of hacky)
ls $in_folder/*HWYTYBGX5*fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 9- | sed -r 's/_R[12].fastq.gz//' \
| sort | uniq | rev | cut -d'_' -f 3- | rev >> $out_folder/base_names.txt

# Now get path to round1 files without R[12].fastq.gz extension
ls $in_folder/*HWYTYBGX5*fastq.gz | sed -r 's/_R[12].fastq.gz//' | uniq | sort >> $out_folder/file_names.txt

# Now get base names for round2 files (flow cell HMYYHBGX7)
ls $in_folder/*HMYYHBGX7*fastq.gz | rev | cut -d'/' -f 1 | rev | cut -d'_' -f 7- | sed -r 's/_R[12].fastq.gz//' \
| sort | uniq | rev | cut -d'_' -f 3- | rev | awk '{print $1"_round2"}' >> $out_folder/base_names.txt

#And finally path to round2 files
ls $in_folder/*HMYYHBGX7*fastq.gz | sed -r 's/_R[12].fastq.gz//' | uniq | sort >> $out_folder/file_names.txt

##### Run in parallel each pair of R1 and R2 reads

# bwa wrapper function
run_bwa(){
bwa mem -t 4 \
-R '@RG\tID:HWYTYBGX5_1\tSM:'"$3"'\tLB:1\tPL:ILLUMINA' \
/workdir/gmv23/bsa/cpepo_genome/Cpepp_v4.1.chr.fa.gz \
${2}_R1.fastq.gz ${2}_R2.fastq.gz | \
samtools view -Sb - > \
${1}/${3}.bam
}
export -f run_bwa

# call wrapper in parallel
parallel --link run_bwa $out_folder \
:::: $out_folder/file_names.txt \
:::: $out_folder/base_names.txt \
