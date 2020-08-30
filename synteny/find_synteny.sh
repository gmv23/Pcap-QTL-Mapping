#/usr/bin/bash

#To make gene names in gff and fasta files line up remove the .1 from fasta files
sed -r 's/>CmoCh(.+)\.1/>CmoCh\1/g' input/Cmoschata_pep_v1.fa > Cmoschata_pep_clean.fa
sed -r 's/>Cp4.1(.+)\.1/>Cp4.1\1/g' input/Cpepo_pep_v4.1.fa > Cpepo_pep_clean.fa

#Make Blast databases
makeblastdb -in Cmoschata_pep_clean.fa -dbtype prot -out Cmoschata
makeblastdb -in Cpepo_pep_clean.fa -dbtype prot -out Cpepo

#Do reciprocal Blast searches
blastp -query Cmoschata_pep_clean.fa -db Cpepo -evalue 1e-5 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out moschata_homs
blastp -query Cpepo_pep_clean.fa -db Cmoschata -evalue 1e-5 -max_target_seqs 5 -max_hsps 1 -outfmt 6 -out pepo_homs

#Make joint gff file
awk '$3 == "gene"' input/Cpepo_gff_v4.1 > Cpepo_genes.gff
paste <(cut -f1 Cpepo_genes.gff | sed -r 's/Cp4.1LG0{,1}(.*)/Cp\1/') \
<(sed -r 's/.*ID=(.*);Name.*/\1/' Cpepo_genes.gff) \
<(cut -f4,5 Cpepo_genes.gff) > Cpepo.bed

awk '$3 == "gene"' input/Cmoschata_gff3_v1 > Cmoschata_genes.gff
paste <(cut -f1 Cmoschata_genes.gff | sed -r 's/Cmo_Chr0{,1}(.*)/Cm\1/') \
<(sed -r 's/.*ID=(.*);Name.*/\1/' Cmoschata_genes.gff) \
<(cut -f4,5 Cmoschata_genes.gff) > Cmoschata.bed

cat Cpepo.bed Cmoschata.bed > squash.gff

#Make joint blast results file
cat pepo_homs moschata_homs > squash.blast

#Put blast and gff files into one directory and run MCScanX
mkdir MC_input
mv squash* MC_input
/workdir/gmv23/software/MCScanX/MCScanX MC_input/squash
