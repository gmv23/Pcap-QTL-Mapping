#!/usr/bin/bash

#This script is to
#Find genes in QTL search regions
#And cross reference with other sources of information on genes to narrow down list of candidates

#Use get_recip_hom.py script to BLAST C. pepo proteins against Cucumis melo proteins
#And find reciprocal best homologs
python /workdir/gmv23/squashQTL/Pcap-QTL-Mapping/CandidateGenes/get_recip_hom.py Cpepo_pep_v4.1.fa CM3.6.1_pep.fasta

#Filter VCF file with SNPs and INDELs to get high quality variants that are segregating
#Also remove variants with MAF<0.20 to get rid of spurious variants that share position with retained position
#MAF of 0.20 means minor allele called in at least 5 samples (pools or parents)
vcftools --vcf ../bsa_bioinformatics_2020-04-30/variants/final_geno.vcf --positions hq_polymorphic_positions.txt --maf 0.20 --recode --out tmp/geno_poly

#Annotate variants
java -jar /workdir/gmv23/bsa/candidates/snpEff-4.3/snpEff.jar Cp4.1 -v \
-c /workdir/gmv23/bsa/candidates/snpEff-4.3/snpEff.config tmp/geno_poly.recode.vcf > tmp/geno_anno.vcf

#Pull out genes from gff and sort
awk '$3 == "gene"' Cpepo_gff_v4.1 | sed 's/\s/\t/g' | sed -r 's/^(.*?)ID=.*?Name=(.*?);/\1\2/' > tmp/gff_genes.gff
bedtools sort -i tmp/gff_genes.gff > tmp/gff_genes_sorted.gff

#Now get list of genes in each region using bedtools
bedtools map -a QTL_Search_Regions.bed -b tmp/gff_genes_sorted.gff -c 9 -o distinct > tmp/genes_in_qtl.txt

#Run python script to print each gene on a different line
python2 /workdir/gmv23/squashQTL/Pcap-QTL-Mapping/CandidateGenes/genes_fat_to_skinny.py tmp/genes_in_qtl.txt tmp/candidate_genes0.txt

#Add gene coordinates from GFF
cut -f 4,5,7,9 tmp/gff_genes_sorted.gff | sort -k 4 > tmp/gff_genes_sorted_by_gene.txt
join tmp/candidate_genes0.txt tmp/gff_genes_sorted_by_gene.txt -1 2 -2 4 > tmp/candidate_genes1.txt

#Look at list of genes with variants as output by snp_eff
#And add counts of moderate and high impact variants in each gene
tail -n +3 snpEff_genes.txt | awk -v OFS='\t' '{print $1, $7, $5}' | sort -k 1 > tmp/variant_counts.txt

#Add variant count information to gene summary
join -e FALSE -a 1 -o 1.1 1.2 1.3 1.4 1.5 2.2 2.3 tmp/candidate_genes1.txt tmp/variant_counts.txt > tmp/candidate_genes2.txt

#Run python script that will cross-reference RBH file and melon DE file to add information 
#on melon homologs and if theyre differentially expressed
python2 /workdir/gmv23/squashQTL/Pcap-QTL-Mapping/CandidateGenes/DE_homolog.py tmp/candidate_genes2.txt tmp/candidate_genes3.txt

#Add gene descriptions
sort -k 1 cucurbita_pepo_gene_description.txt | sed 's/\s/\t/' > tmp/sorted_gene_descriptions.txt
join -t $'\t' -e NA -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.2 tmp/candidate_genes3.txt tmp/sorted_gene_descriptions.txt > tmp/candidate_genes4.txt

#Add KEGG annotations
sort -k 1 Cpepo_v4.1_KO.txt > tmp/sorted_kegg.txt
join -t $'\t' -e NA -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.3 tmp/candidate_genes4.txt tmp/sorted_kegg.txt > tmp/candidate_genes5.txt

#Join GO terms for each gene with ; and then add GO terms
python /workdir/gmv23/squashQTL/Pcap-QTL-Mapping/CandidateGenes/go_skinny_to_fat.py Cpepo_GO_anno.txt tmp/go_fat.txt
sort -k 1 tmp/go_fat.txt > tmp/go_fat_sorted.txt
join -t $'\t' -e NA -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.3 tmp/candidate_genes5.txt tmp/go_fat_sorted.txt > candidate_genes_list.txt
