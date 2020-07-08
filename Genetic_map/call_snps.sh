#!/usr/bin/bash

#This script is to call SNPs from paired end GBS reads of F2 squash plants using GBS-SNP-CROP

#Steps 1-3 must be run separately for each library
#Step 1 must be run twice to recognize both CAGC and CTGC ApeK1 residues

#Path to directory containg GBS-SNP-CROP scripts
path_to_gsc="/workdir/gmv23/software/GBS-SNP-CROP/GBS-SNP-CROP-scripts/v.4.1"

########### Step 1: Parsing reads, trimming barcodes cut sites and adapters, appending barcodes to read headers ###############

#Plates 1 and 2 in parallel, cut site CAGC
parallel --link perl $path_to_gsc/GBS-SNP-CROP-1.pl -d PE -b {1} \
-fq {2} -s 1 -e 1 -enz1 CAGC -enz2 CAGC -t 24 :::: barcode_paths.txt :::: fastqfileseeds.txt

#Move all output files to separate directories so they will not be rewritten for cut site CTGC
mkdir cutsite_CAGC
mv distribs/ cutsite_CAGC/
mv parsed/ cutsite_CAGC/
mv singles/ cutsite_CAGC/ 
mv summaries/ cutsite_CAGC/

#Now do Plates 1 and 2 in parallel, for cut site CTGC
parallel --link perl $path_to_gsc/GBS-SNP-CROP-1.pl -d PE -b {1} \
-fq {2} -s 1 -e 1 -enz1 CTGC -enz2 CTGC -t 24 :::: barcode_paths.txt :::: fastqfileseeds.txt

#Move these output files to cut site CTGC specific directories
mkdir cutsite_CTGC
mv distribs/ cutsite_CTGC/
mv parsed/ cutsite_CTGC/
mv singles/ cutsite_CTGC/ 
mv summaries/ cutsite_CTGC/

#Now concatenate CTGC reads and CAGC reads into same file
mkdir parsed
cat <(gunzip -c cutsite_CAGC/parsed/17G53-Plate1_S32_L001_001.R1parsed.fq.gz) <(gunzip -c cutsite_CTGC/parsed/17G53-Plate1_S32_L001_001.R1parsed.fq.gz) | gzip > parsed/17G53-Plate1_S32_L001_001.R1parsed.fq.gz
cat <(gunzip -c cutsite_CAGC/parsed/17G53-Plate1_S32_L001_001.R2parsed.fq.gz) <(gunzip -c cutsite_CTGC/parsed/17G53-Plate1_S32_L001_001.R2parsed.fq.gz) | gzip > parsed/17G53-Plate1_S32_L001_001.R2parsed.fq.gz
cat <(gunzip -c cutsite_CAGC/parsed/17G53-Plate2_S33_L001_001.R1parsed.fq.gz) <(gunzip -c cutsite_CTGC/parsed/17G53-Plate2_S33_L001_001.R1parsed.fq.gz) | gzip > parsed/17G53-Plate2_S33_L001_001.R1parsed.fq.gz
cat <(gunzip -c cutsite_CAGC/parsed/17G53-Plate2_S33_L001_001.R2parsed.fq.gz) <(gunzip -c cutsite_CTGC/parsed/17G53-Plate2_S33_L001_001.R2parsed.fq.gz) | gzip > parsed/17G53-Plate2_S33_L001_001.R2parsed.fq.gz

###################################################  Step 2: Trimmomatic wrapper  ######################################################

#Note: Comment out line in GBS-SNP-CROP-2.pl that deletes all parsed read files prior to trimming once completed; otherwise any problems at this stage in script will result in deletion of previous step output

cd parsed

#Copy illumina adaptor file to directory
cp /programs/trimmomatic/adapters/TruSeq3-PE.fa .

#Run trimmomatic in parallel on plates 1 and 2 using default parameters supplied in GBS-SNP-CROP-2.pl
parallel perl $path_to_gsc/GBS-SNP-CROP-2.pl -tm /programs/trimmomatic/trimmomatic-0.39.jar -d PE -fq {} -t 24 -ad TruSeq3-PE.fa:2:30:10:8:true :::: ../fastqfileseeds.txt

#########################################################  Step 3: Demultiplex  ######################################################

#Note: Comment out lines in GBS-SNP-CROP-3.pl that create demultiplexed directory and raise an exception if directory already made
#Otherwise you will not be able to run on second plate since directory already exists from launching script for first plate
mkdir demultiplexed

#Run step 3 in parallel on plates 1 and 2
parallel --link perl $path_to_gsc/GBS-SNP-CROP-3.pl -d PE -b ../{1} -fq {2} :::: ../barcode_paths.txt :::: ../fastqfileseeds.txt

######################################  Step 5: Align reads to reference and make mpileup files  ###########################################

#Note: Step 4 skipped because we are using reference genome

#Change directory and bring reference genome to directory
cd demultiplexed
cp /workdir/gmv23/bsa/cpepo_genome/Cpepp_v4.1.chr.fa .

#Use default alignment filtering parameters provided in GBS-SNP-CROP-5.pl
#Use combined_barcodeID.txt file which has placefiller for barcodes -- at this point we are using samples from two different libraries that have same barcode
perl $path_to_gsc/GBS-SNP-CROP-5.pl -bw /programs/bwa-0.7.17/bwa -st /programs/samtools-1.9/bin/samtools \
-d PE -b ../../barcodes/combined_barcodeID.txt -ref Cpepp_v4.1.chr.fa -t 36

#################################################  Step 6: Create variant master matrix  ##################################################
perl $path_to_gsc/GBS-SNP-CROP-6.pl -b ../../barcodes/combined_barcodeID.txt -out SquashF2.MasterMatrix.txt -t 36

#################################################  Step 7: Filter variants and call genotypes  ##################################################
perl $path_to_gsc/GBS-SNP-CROP-7.pl -in SquashF2.MasterMatrix.txt -out SquashF2.GenoMatrix.txt \
-mnHoDepth0 5 -mnHoDepth1 20 -mnHetDepth 3 -altStrength 0.9 -mnAlleleRatio 0.15 -mnCall 0.75 -mnAvgDepth 7 -mxAvgDepth 150
#The parameters that are not default are altStrength, mnAlleleRatio, mnAvgDepth, mxAvg Depth

#################################################  Step 8: Convert to dosage matrix  ##################################################
cd variants
#Convert to 0,0.5,1 dosage matrix 
perl $path_to_gsc/GBS-SNP-CROP-8.pl -in SquashF2.GenoMatrix.txt -out SquashF2 -b ../../../barcodes/combined_barcodeID.txt -formats R

#Add header row and some extra header columns for filtering purposes
inds=$(cut -f2 ../../../barcodes/combined_barcodeID.txt | tr '\n' '\t')
head="CHROM\tBP\tTYPE\tREF\tALT\tREAD_DEPTH\tCALL_RATE\tHOMO_REF\tHET\tHOMO_ALT\tMAF\tZSCORE\t"
header="$head""$inds"
paste <(cut -f1-12 SquashF2.GenoMatrix.txt) <(cut -f3- SquashF2.R.txt) > tmp.txt
echo -e "$header" | cat - tmp.txt > SquashF2_geno.txt
rm tmp.txt
