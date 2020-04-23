#!/usr/bin/bash

#Align reads to genome in parallel

reads_folder=$1

# Make files with base names, file paths, and flow cell IDs  for input to bwa in parallel

# Get base names from round1 files (flow cell HWYTYBGX5) (Kind of hacky)
ls $reads_folder/*HWYTYBGX5*fastq.gz | rev | cut -d'/' -f 1 | rev | \
cut -d'_' -f 9- | sed -r 's/_R[12].fastq.gz//' | sort | uniq | \
rev | cut -d'_' -f 3- | rev >> base_names.txt

# Now get path to round1 files without R[12].fastq.gz extension
ls $reads_folder/*HWYTYBGX5*fastq.gz | sed -r 's/_R[12].fastq.gz//' | \
sort | uniq >> file_names.txt

# Now get flowcell for round1 files
ls $reads_folder/*HWYTYBGX5*fastq.gz | sed -r 's/_R[12].fastq.gz//' | \
sort | uniq | rev | cut -d '/' -f 1 | rev | \
cut -d '_' -f 4 >> flowcell_names.txt

# Now get base names for round2 files (flow cell HMYYHBGX7)
ls $reads_folder/*HMYYHBGX7*fastq.gz | rev | cut -d'/' -f 1 | rev | \
cut -d'_' -f 7- | sed -r 's/_R[12].fastq.gz//' | sort | uniq | \
rev | cut -d'_' -f 3- | rev | awk '{print $1"_round2"}' >> base_names.txt

#And path to round2 files
ls $reads_folder/*HMYYHBGX7*fastq.gz | sed -r 's/_R[12].fastq.gz//' | \
sort | uniq >> file_names.txt

#And flowcell for round2 files
ls $reads_folder/*HMYYHBGX7*fastq.gz | sed -r 's/R[12].fastq.gz//' | \
sort | uniq | rev | cut -d '/' -f 1 | rev | \
cut -d '_' -f 4 >> flowcell_names.txt

##### Run in parallel each pair of R1 and R2 reads

# bwa wrapper function
run_bwa(){
bwa mem -t 4 \
-R '@RG\tID:'"$3"'_1\tSM:'"$2"'\tLB:1\tPL:ILLUMINA' \
/workdir/gmv23/bsa/cpepo_genome/Cpepp_v4.1.chr.fa.gz \
${1}_R1.fastq.gz ${1}_R2.fastq.gz | \
samtools view -Sb - > \
bam_files/${2}.bam
}
export -f run_bwa

# call wrapper in parallel
parallel --link --jobs 20 run_bwa \
:::: file_names.txt \
:::: base_names.txt \
:::: flowcell_names.txt
