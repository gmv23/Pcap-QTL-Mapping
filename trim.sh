#!/usr/bin/bash

#Run fastp in parallel

folder=$1
out_folder=$2

if [ ! -d $2/qc ]; then
	mkdir $2/qc
fi

ls $folder/*fastq.gz | sed -r 's/_R[12].*//' | uniq | \
parallel /workdir/gmv23/bsa/fastp_current/fastp -i {}_R1.fastq.gz -I {}_R2.fastq.gz  \
-h $out_folder/qc/{/.}.html \
-q 15 -u 50 \
-l 20 \
--cut_by_quality5 --cut_by_quality3 -W 4 -M 20 \
--correction \
--trim_poly_g \
-w 3 \
-o $out_folder/{/}_trimmed_R1.fastq.gz \
-O $out_folder/{/}_trimmed_R2.fastq.gz

