setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/BSA/")

######################################################### Import data and load scripts  #############################################################

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
convert_alleles <- function(s_geno, t_geno, bulk_alleles){
  s_counts <- rep(NA, length(s_geno))
  t_counts <- rep(NA, length(t_geno))
  s_counts[s_geno < t_geno] <- bulk_alleles[s_geno < t_geno, 1]
  s_counts[s_geno > t_geno] <- bulk_alleles[s_geno > t_geno, 2]
  t_counts[s_geno < t_geno] <- bulk_alleles[s_geno < t_geno, 2]
  t_counts[s_geno > t_geno] <- bulk_alleles[s_geno > t_geno, 1]
  return(cbind(s_counts,t_counts))
}

for(i in seq(3,ncol(counts.bi),by=2)){
  counts.bi[,i:(i+1)] <- convert_alleles(ps_geno, pt_geno, counts.bi[,i:(i+1)])
}

####################################################### Compare tech reps and collapse #############################################################

#Kind of complicated number of libraries and lanes per each pool
#Good way to compare variation due to rep vs library vs sequencer?

#For now just collapse all the tech reps, throw away RAN-1 and TOL-2 merge
counts.bi.collapse <- data.frame("CHROM" = counts.bi$CHROM,
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

#####################################   Function to filter within pools and convert counts to frequencies   #################################

#Filter counts and get allele frequencies or read depth
#Used in collapse_columns script
filter_and_calculate <- function(counts, rd_thresh, freq_thresh, stat_type = c("frequency", "read depth")){
  total <- sum(counts)
  if(total == 0){
    return(NA)
  }
  else{
    af <- counts[2]/total
    if(total < rd_thresh[1] | total > rd_thresh[2] | af > freq_thresh | af < (1-freq_thresh)){
      return(NA)
    }else{
      if(stat_type == "frequency"){
        return(counts[2]/total)
      }else if(stat_type == "read depth"){
        return(sum(counts))
      }else{
        stop("Invalid stat_type")
      }
    }
  }
}

#Turn A1 and A2 count columns into one column ---- either total read depth or frequency
collapse_columns <- function(counts, data_columns, 
                             rd_filter, rd_filter_type = c("absolute", "quantile"),
                             freq_filter, stat_type = c("frequency", "read depth")){
  
  out <- data.frame("CHROM" = counts[,1],
                    "POS" = counts[,2],
                    matrix(NA,nrow=nrow(counts), ncol=length(data_columns)/2))
  sample_starts <- data_columns[seq(1,length(data_columns),by=2)]
  fill_cols <- 3:ncol(out)
  
  for(i in 1:length(fill_cols)){
    
    #Get minimum read depth for each pool
    pool_counts <- counts[,sample_starts[i]:(sample_starts[i]+1)]
    cov_dist <- apply(pool_counts,1, sum)
    if(rd_filter_type == "quantile"){
      rd_thresh <- quantile(cov_dist, rd_filter)
    }else if(rd_filter_type == "absolute"){
      rd_thresh <- rd_filter
    }else{
      stop("rd_filter must be either 'absolute' or quantile' ")
    }
    out[,fill_cols[i]] <- apply(pool_counts, 1, filter_and_calculate,
                                rd_thresh = rd_thresh,
                                freq_thresh = freq_filter, 
                                stat_type = stat_type)
  }
  
  #Rename colnames assuming they are in format "pool_A[12]"
  colnames.base <- sapply(colnames(counts)[data_columns], function(x) unlist(substr(x, 1, nchar(x)-3)))
  colnames(out)[fill_cols] <- colnames.base[seq(1,length(colnames.base),2)]
  return(out)
}
############################################################ Export files  #############################################################

allele_freqs <- collapse_columns(counts = counts.bi.collapse, data_columns = 3:ncol(counts.bi.collapse), stat_type = "frequency",
                                 rd_filter = c(0.10,1), rd_filter_type = "quantile", freq_filter = 0.9)

# Write marker files for multipool -- separate file for each pool and each chromosome
# Markers will be filtered within pools and markers filtered within that pool will not be included in any comparison with that pool
# Therefore there will be 8 files produced for each chromosome --- one for each pool for each of 4 comparisons

for(chrom in unique(counts.bi.collapse$CHROM)){
  chrom <- as.character(chrom)
  chrom_num <- unlist(strsplit(chrom, "Cp4.1LG"))[2]
  
  #S1 v T1
  rows_to_remove <- apply(allele_freqs[,c('s1', 't1')],1, function(x) any(is.na(x)))
  counts.filter <- counts.bi.collapse[!rows_to_remove,]
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "s1_A1", "s1_A2")], 
              paste('data/multipool/s1_v_t1/chr', chrom_num, '_S1.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "t1_A1", "t1_A2")], 
              paste('data/multipool/s1_v_t1/chr', chrom_num, '_T1.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  
  #S2 v T2
  rows_to_remove <- apply(allele_freqs[,c('s2', 't2')],1, function(x) any(is.na(x)))
  counts.filter <- counts.bi.collapse[!rows_to_remove,]
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "s2_A1", "s2_A2")], 
              paste('data/multipool/s2_v_t2/chr', chrom_num, '_S2.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "t2_A1", "t2_A2")], 
              paste('data/multipool/s2_v_t2/chr', chrom_num, '_T2.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  
  #S1 v S2
  rows_to_remove <- apply(allele_freqs[,c('s1', 's2')],1, function(x) any(is.na(x)))
  counts.filter <- counts.bi.collapse[!rows_to_remove,]
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "s1_A1", "s1_A2")], 
              paste('data/multipool/s1_v_s2/chr', chrom_num, '_S1.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "s2_A1", "s2_A2")], 
              paste('data/multipool/s1_v_s2/chr', chrom_num, '_S2.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  
  #T1 v T2
  rows_to_remove <- apply(allele_freqs[,c('t1', 't2')],1, function(x) any(is.na(x)))
  counts.filter <- counts.bi.collapse[!rows_to_remove,]
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "t1_A1", "t1_A2")], 
              paste('data/multipool/t1_v_t2/chr', chrom_num, '_T1.txt', sep=""),
              quote = F, row.names = F, col.names = F)
  write.table(counts.filter[counts.filter$CHROM==chrom,c("POS", "t2_A1", "t2_A2")], 
              paste('data/multipool/t1_v_t2/chr', chrom_num, '_T2.txt', sep=""),
              quote = F, row.names = F, col.names = F)
}

#Now write file with filtered allele frequencies
write.table(allele_freqs, "data/freqs_filtered.txt", quote=F, row.names = F, col.names = T)
