setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/BSA/")

######################################################### Import data and load scripts  #############################################################

counts <- read.table("data/allele_counts.txt", stringsAsFactors = F, header=T)
nrow(counts)
seq_stats <- read.csv("data/sequencing_stats.txt", stringsAsFactors=F, header=F)
colnames(seq_stats) <- c("Sample", "Raw", "Filtered", "Aligned")
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

####################################################### Compare tech reps and collapse #############################################################

#Kind of complicated number of libraries and lanes per each pool

#Turn allele counts into allele freqs for all tech reps separately
allele_freqs.full <- collapse_columns(counts = counts.bi, data_columns = 3:ncol(counts.bi), stat_type = "frequency",
                         rd_filter = c(0.10,1), rd_filter_type = "quantile", freq_filter = 0.9)

#Get pairwise correlations of allele frequencies
af.cor <- cor(allele_freqs.full[-c(1,2)],use="complete.obs")
af.cor[lower.tri(af.cor)] <- NA
af.cor[diag(nrow=nrow(af.cor))==1] <- NA

#Generate list of sample pairs representing either tech reps or bio reps
tech_reps <- rbind(c(3,4),
                   c(6,7),
                   c(9,10),
                   c(11,12),
                   c(11,13),
                   c(11,14),
                   c(12,13),
                   c(13,14),
                   c(12,14),
                   c(15,16))
bio_reps <- rbind(c(3,11),
                  c(3,12),
                  c(3,13),
                  c(3,14),
                  c(4,11),
                  c(4,12),
                  c(4,13),
                  c(4,14),
                  c(6,15),
                  c(6,16),
                  c(7,15),
                  c(7,16))

#Pull out read depth medians of all pools
rds <- collapse_columns(counts = counts.bi, data_columns = 3:ncol(counts.bi), stat_type = "read depth",
                              rd_filter = c(0.10,1), rd_filter_type = "quantile", freq_filter = 0.9)
rds_medians <- apply(rds[3:ncol(rds)], 2, median, na.rm=T)

#Now get data frame showing pairwise correlations between all pairs of samples
#and include whether those samples are a bio or tech rep or neither
#and what the avg median read depth is
frequency_correlations <- data.frame("Sample1" = rep(NA,16*16/2-16),
                                     "Sample2" = NA,
                                     "Relation" = NA,
                                     "Cor" = NA,
                                     "RD1" = NA,
                                     "RD2" = NA,
                                     "RDmean" = NA)

counter <- 0
for(i in 1:15){
  for(j in (i+1):16){
    counter <- counter + 1
    sample1 <- colnames(af.cor)[i]
    sample2 <- colnames(af.cor)[j]
    samplecor <- af.cor[i,j]
    rd1 <- rds_medians[i]
    rd2 <- rds_medians[j]
    RDmean <- mean(c(rd1,rd2))
    frequency_correlations[counter,] <- c(sample1, sample2, "NA", samplecor, rd1, rd2, RDmean)
  }
}

for(i in 1:nrow(tech_reps)){
  frequency_correlations$Relation[frequency_correlations$Sample1==colnames(af.cor)[tech_reps[i,1]] & 
                                    frequency_correlations$Sample2==colnames(af.cor)[tech_reps[i,2]]] <- "TECH"
}
for(i in 1:nrow(bio_reps)){
  frequency_correlations$Relation[frequency_correlations$Sample1==colnames(af.cor)[bio_reps[i,1]] & 
                                    frequency_correlations$Sample2==colnames(af.cor)[bio_reps[i,2]]] <- "BIO"
}

#Plot how average pool read depth affects correlation between allele frequency values between pools
#color data points by whether they're a bio rep, tech rep, or neither

frequency_correlations$Cor <- as.numeric(frequency_correlations$Cor)
frequency_correlations$RD1 <- as.numeric(frequency_correlations$RD1)
frequency_correlations$RD2 <- as.numeric(frequency_correlations$RD2)
frequency_correlations$RDmean <- as.numeric(frequency_correlations$RDmean)
relation_color <- rep('gray',nrow(frequency_correlations))
relation_color[frequency_correlations$Relation=="TECH"] <- 'blue'
relation_color[frequency_correlations$Relation=="BIO"] <- 'purple'
plot(frequency_correlations$RDmean, frequency_correlations$Cor, col=relation_color, pch=16)
cor(frequency_correlations$RDmean, frequency_correlations$Cor)
cor(frequency_correlations$RDmean[frequency_correlations$Relation=="TECH"], 
    frequency_correlations$Cor[frequency_correlations$Relation=="TECH"])
cor(frequency_correlations$RDmean[frequency_correlations$Relation=="BIO"], 
    frequency_correlations$Cor[frequency_correlations$Relation=="BIO"])
mean(frequency_correlations$Cor[frequency_correlations$Relation=="TECH"])
mean(frequency_correlations$Cor[frequency_correlations$Relation=="BIO"])

#Now just collapse all the tech reps, throw away RAN-1 and TOL-2 merge
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

counts.parents <- data.frame("CHROM" = counts.bi$CHROM,
                             "POS" = counts.bi$POS,
                             "D25_A1" = counts.bi$B02_D25_A1,
                             "D25_A2" = counts.bi$B02_D25_A2,
                             "156015_A1" = counts.bi$C02_15_6015_A1,
                             "156015_A2" = counts.bi$C02_15_6015_A2)


###################################################### Make sequencing stats table  #############################################################

rds_pools <- collapse_columns(counts = counts.bi.collapse, data_columns = 3:ncol(counts.bi.collapse), stat_type = "read depth",
                              rd_filter = c(0.10,1), rd_filter_type = "quantile", freq_filter = 0.9)
rds_parents <- collapse_columns(counts = counts.parents, data_columns = 3:ncol(counts.parents), stat_type = "read depth",
                                rd_filter = c(0,1), rd_filter_type = "quantile", freq_filter = 1)

rds_pools_medians <- apply(rds_pools[3:ncol(rds_pools)], 2, median, na.rm=T)
rds_parents_medians <- apply(rds_parents[3:ncol(rds_parents)], 2, median)

genomesize <- 263500453
Gb <- 1000000000

seq_table <- data.frame("Sample" = seq_stats$Sample,
                        "Raw reads (Gb)" = round(seq_stats$Raw/Gb,2),
                        "Raw reads (Genome coverage)" = round(seq_stats$Raw/genomesize,0),
                        "Filtered reads (Gb)" = round(seq_stats$Filtered/Gb,2),
                        "Filtered reads (Genome coverage)" = round(seq_stats$Filtered/genomesize,0),
                        "Aligned reads (Gb)" = round(seq_stats$Aligned/Gb,2),
                        "Aligned reads (Genome coverage)" = round(seq_stats$Aligned/genomesize,0),
                        "Median SNP read depth" = c(rds_parents_medians, rds_pools_medians)
)

write.csv(seq_table, "tables/seqstats.csv", row.names = F, quote=F)


############################################################ Export files  #############################################################

allele_freqs <- collapse_columns(counts = counts.bi.collapse, data_columns = 3:ncol(counts.bi.collapse), stat_type = "frequency",
                                 rd_filter = c(0.10,1), rd_filter_type = "quantile", freq_filter = 0.9)

n_snps <- apply(allele_freqs[,3:ncol(allele_freqs)],2,function(x) sum(!is.na(x)))
summary(n_snps)

# Write marker files for multipool -- separate file for each pool and each chromosome
# Markers will be filtered within pools and markers filtered within that pool will not be included in any comparison with that pool
# Therefore there will be 8 files produced for each chromosome --- one for each pool for each of 4 comparisons

#Save table with 'unfiltered' allele counts (unfiltered within pools, but filtering steps that apply to all pools already performed)
counts.clean <- counts.bi.collapse
colnames(counts.clean) <- c("CHROM", "POS", "SUS1_A1", "SUS1_A2", "RES1_A1", "RES1_A2", "SUS2_A1", "SUS2_A2", "RES2_A1", "RES2_A2", "RAN2_A1", "RAN2_A2")
write.table(counts.clean, 
            'data/filtered_allele_counts.txt',
            quote = F, row.names = F, col.names = T, sep="\t")

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
