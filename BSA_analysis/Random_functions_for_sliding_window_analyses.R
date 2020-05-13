##################  Convert ref and alternative to susceptible-derived and tolerant-derived   ###################

convert_alleles <- function(s_geno, t_geno, bulk_alleles){
  s_counts <- rep(NA, length(s_geno))
  t_counts <- rep(NA, length(t_geno))
  s_counts[s_geno < t_geno] <- bulk_alleles[s_geno < t_geno, 1]
  s_counts[s_geno > t_geno] <- bulk_alleles[s_geno > t_geno, 2]
  t_counts[s_geno < t_geno] <- bulk_alleles[s_geno < t_geno, 2]
  t_counts[s_geno > t_geno] <- bulk_alleles[s_geno > t_geno, 1]
  return(cbind(s_counts,t_counts))
}

###############################   Filter counts and convert to frequencies   #################################

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

#################################  Find pools where depth of A1 or A2 is 0  ##############################

find_zero_depth <- function(counts, data_columns){
  zero_depth_rows <- c()
  for(i in data_columns[seq(1,length(data_columns), by=2)]){
     zero_depth_rows <- c(zero_depth_rows,
                          which(apply(counts[,i:(i+1)], 1, function(x) all(x==0))))
  }
  return(unique(zero_depth_rows))
}

##############################  Simulation of confidence intervals for deltaSNP test  ##########################

simulate_confidence <- function(rds, perms=1000, n, alpha=0.05){
  
  rd1 <- rds[1]
  rd2 <- rds[2]
  
  #First simulate sampling of segregants in bulks
  #binomial sampling of 2 alleles per individual, allele frequency at 0.5
  pb1 <- rbinom(prob = 0.5, n = perms, size = n*2)/(n*2)
  pb2 <- rbinom(prob = 0.5, n = perms, size = n*2)/(n*2)
  
  #Now use binomial sampling with read depth to get snp indices
  snpin1 <- rbinom(prob=pb1, n = perms, size = rd1)/rd1
  snpin2 <- rbinom(prob=pb2, n = perms, size = rd2)/rd2
  deltasnp <- snpin1 - snpin2
  
  lower <- quantile(deltasnp, probs = alpha/2, na.rm=T)
  upper <- quantile(deltasnp, probs = 1-alpha/2, na.rm = T)
  
  return(c(lower,upper))
}



