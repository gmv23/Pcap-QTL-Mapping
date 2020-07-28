setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/CandidateGenes/")

#This script is modified from Filter_sites.R
#It loads a vcf file with both SNPs and INDELS on BSA pools
#And finds quality sites that are either fixed for opposite alleles in the parents
#Or homo in one parent and het in the other
#This is to look for possible deleterious variants segregating in candidate genes

######################################################### Import data  #############################################################
counts <- read.table("allele_counts_all_variants.txt", stringsAsFactors = F, header=T)
nrow(counts)
############################################################ Filter and recode data #############################################################

#Get rid of chromosome 00
counts <- counts[counts$CHROM!="Cp4.1LG00",]

######## Filter SNPs based on quality -- read depth and MAF

#Find spurious SPNs -----> how many have very low counts of minor allele
A1_columns <- grep("A1", colnames(counts))
A2_columns <- grep("A2", colnames(counts))
# minor allele count frequency
A1_counts <- apply(counts[,A1_columns], 1, sum)
A2_counts <- apply(counts[,A2_columns], 1, sum)
macf <- apply(cbind(A1_counts, A2_counts), 1, function(x)
  min(x)/sum(x))
total_reads <- apply(cbind(A1_counts, A2_counts), 1, sum)
hist(macf)
hist(total_reads)
#Filter on: total read depth < 95th percentile and > 5th percentile; MAF > 0.1
macf_threshold <- 0.1
total_rd_max <- quantile(total_reads, 0.95)
total_rd_min <- quantile(total_reads, 0.05)
# 'high quality' counts after filtering
counts.hq <- counts[macf > macf_threshold & 
                      total_reads < total_rd_max &
                      total_reads > total_rd_min,]
nrow(counts.hq)

### Now find SNPs that are either homozygous for opposite alleles in parents
### or heterozygous in one parent and homozygous in the other
counts.poly <- counts.hq[( (counts.hq$B02_D25_geno==2 & counts.hq$C02_15_6015_geno==0) |
                            (counts.hq$B02_D25_geno==0 & counts.hq$C02_15_6015_geno==2) |
                            (counts.hq$B02_D25_geno==1 & (counts.hq$C02_15_6015_geno==0 | counts.hq$C02_15_6015_geno==2)) |
                            ((counts.hq$B02_D25_geno==0 | counts.hq$B02_D25_geno==2) & counts.hq$C02_15_6015_geno==1)),1:2]
nrow(counts.poly)
#Write table of chromosomes and positions
write.table(counts.poly, "hq_polymorphic_positions.txt", row.names=F, quote=F, col.names=F)

