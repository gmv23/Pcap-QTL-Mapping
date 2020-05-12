#!/usr/bin/bash

#Run fastp in parallel

in_file_names=$1
out_folder=$2


cat $in_file_names | sed -r 's/_R[12].*//' | uniq | \
parallel /programs/fastp-0.20.0/bin/fastp -i {}_R1.fastq.gz -I {}_R2.fastq.gz  \
-h $out_folder/qc/{/.}.html \
-q 15 -u 40 \
-l 15 \
--cut_front --cut_tail -W 4 -M 20 \
--correction \
--trim_poly_g \
-w 3 \
-o $out_folder/{/}_trimmed_R1.fastq.gz \
-O $out_folder/{/}_trimmed_R2.fastq.gz

