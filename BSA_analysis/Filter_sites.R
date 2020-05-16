setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/BSA/")

######################################################### Import data and load scripts  #############################################################

source("../scripts/Pcap-QTL-Mapping/BSA_analysis/Random_functions_for_sliding_window_analyses.R")

counts <- read.table("data/allele_counts.txt", stringsAsFactors = F, header=T)
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

### Now filter for SNPs that are fixed for opposite alleles in parents so we can figure out where alleles are derived from
# Just sites where parental read depth > 6 and parents are fixed for opposite alleles
counts.bi <- counts.hq[counts.hq$B02_D25_geno!=1 & counts.hq$C02_15_6015_geno!=1 & 
                         counts.hq$B02_D25_geno!=counts.hq$C02_15_6015_geno &
                         counts.hq$B02_D25_A1 + counts.hq$B02_D25_A2 > 6 
                       & counts.hq$C02_15_6015_A1 + counts.hq$C02_15_6015_A2 > 6,]

#Pull out parental genotypes and remove them and allele counts from data frame
ps_geno <- counts.bi$B02_D25_geno
pt_geno <- counts.bi$C02_15_6015_geno
counts.bi$B02_D25_geno <- NULL
counts.bi$C02_15_6015_geno <- NULL
counts.bi$C02_15_6015_A1 <- NULL
counts.bi$C02_15_6015_A2 <- NULL
counts.bi$B02_D25_A1 <- NULL
counts.bi$B02_D25_A2 <- NULL

#Change A1 and A2 to susceptible parent derived (A1) and tolerant parent derived (A2)
for(i in seq(3,ncol(counts.bi),by=2)){
  counts.bi[,i:(i+1)] <- convert_alleles(ps_geno, pt_geno, counts.bi[,i:(i+1)])
}

####################################################### Compare tech reps and collapse #############################################################

#Kind of complicated number of libraries and lanes per each pool
#Good way to compare variation due to rep vs library vs sequencer?

#For now just collapse all the tech reps, throw away RAN-1 and TOL-2 merge
counts.bi.c1 <- data.frame("CHROM" = counts.bi$CHROM,
                           "POS" = counts.bi$POS,
                           "s1_A1" = counts.bi$B01_1_S_round2_A1 + counts.bi$B01_1_S_A1,
                           "s1_A2" = counts.bi$B01_1_S_round2_A2 + counts.bi$B01_1_S_A2,
                           "t1_A1" = counts.bi$C01_1_T_round2_A1 + counts.bi$C01_1_T_A1,
                           "t1_A2" = counts.bi$C01_1_T_round2_A2 + counts.bi$C01_1_T_A2,
                           "s2_A1" = counts.bi$F01_2_S_1_round2_A1 + counts.bi$F01_2_S_1_A1 + counts.bi$G01_2_S_2_round2_A1 + counts.bi$G01_2_S_2_A1,
                           "s2_A2" = counts.bi$F01_2_S_1_round2_A2 + counts.bi$F01_2_S_1_A2 + counts.bi$G01_2_S_2_round2_A2 + counts.bi$G01_2_S_2_A2,
                           "t2_A1" = counts.bi$A02_2_T_2_round2_A1 + counts.bi$H01_2_T_1_round2_A1 + counts.bi$H01_2_T_1_A1,
                           "t2_A2" = counts.bi$A02_2_T_2_round2_A2 + counts.bi$H01_2_T_1_round2_A2 + counts.bi$H01_2_T_1_A2,
                           "r2_A1" = counts.bi$D01_2_R_1_A1 + counts.bi$E01_2_R_2_A1,
                           "r2_A2" = counts.bi$D01_2_R_1_A2 + counts.bi$E01_2_R_2_A2
)

############################################################ Export files  #############################################################

# Write marker files for multipool -- separate file for each pool and each chromosme
for(chrom in unique(counts.bi.c1$CHROM)){
  chrom <- as.character(chrom)
  chrom_num <- unlist(strsplit(chrom, "Cp4.1LG"))[2]

  #S1
  write.table(counts.bi.c1[counts.bi.c1$CHROM==chrom,c(2,3,4)], 
              paste('data/multipool/chr', chrom_num, '_S1.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  #T1
  write.table(counts.bi.c1[counts.bi.c1$CHROM==chrom,c(2,5,6)], 
              paste('data/multipool/chr', chrom_num, '_T1.txt', sep=""),
                    quote = F, row.names = F, col.names = F)
  #S2
  write.table(counts.bi.c1[counts.bi.c1$CHROM==chrom,c(2,7,8)],
              paste('data/multipool/chr', chrom_num, '_S2.txt', sep=""),
                    quote = F, row.names = F, col.names = F)
  #T2
  write.table(counts.bi.c1[counts.bi.c1$CHROM==chrom,c(2,9,10)], 
              paste('data/multipool/chr', chrom_num, '_T2.txt', sep=""),
                    quote = F, row.names = F, col.names = F)
}

#Also write filtered, re-coded, and collapsed counts file
write.table(counts.bi.c1, "data/counts_filtered.txt", quote=F, row.names = F, col.names = T)
