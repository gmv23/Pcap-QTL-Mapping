setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/BSA/")
library(reshape2)
library(grid)

#########################################################   Import data and load scripts  #############################################################

allele_freqs <- read.table("data/freqs_filtered.txt", stringsAsFactors = F, header=T)
lods <- read.table("data/multipool_lods.txt", header=T)
scaffold_sizes <- read.table("data/scaffold_sizes.txt")
multi_peaks <- read.csv("data/CIs_and_peaks.txt")

####################################################### Process data prior to making plot ###################################

#Add column names to scaffold sizes and remove LG 00
colnames(scaffold_sizes) <- c("CHROM", "Length")
scaffold_sizes <- scaffold_sizes[scaffold_sizes$CHROM!="Cp4.1LG00",]
scaffold_sizes$CHROM <- droplevels(scaffold_sizes$CHROM)

#Get information out of scaffold sizes matrix
lgs.unique <- scaffold_sizes[,1]
lgs.length <- scaffold_sizes[,2]

#Make cumulative positions required to plot all chromosomes in one plot for LOD score plot
lgs <- lods$CHROM
lgs.cumsum <- c(0,cumsum(scaffold_sizes[,2]))
cumpos <- rep(NA, nrow(lods))
for(i in 1:length(lgs.unique)){
  scaff <- lgs.unique[i]
  add <- lgs.cumsum[i]
  cumpos[lods$CHROM==scaff] <- lods$POS[lods$CHROM==scaff] + add
}

#Make matrix with coordinates of chromosome centers for tick marks
lgs.centers <- matrix(NA, nrow=length(lgs.unique), ncol=2)
for(i in 1:nrow(lgs.centers)){
  lgs.centers[i,1] <- as.character(lgs.unique)[i]
  lgs.centers[i,2] <- (lgs.cumsum[i] + lgs.cumsum[i+1])/2
}

#Make coordinates of 5 MB intervals for tick mark positions 
tick_names <- c()
tick_cumpos <- c()
tick_spacing <- 1000000
for(i in 1:length(lgs.unique)){
  scaff <- lgs.unique[i]
  pos.lg <- lods$POS[lods$CHROM==scaff]
  lg.tick_names <- seq(0,max(pos.lg), by = tick_spacing)[-1]
  lg.tick_cumpos <- lg.tick_names + lgs.cumsum[i]
  tick_names <- c(tick_names, lg.tick_names)
  tick_cumpos <- c(tick_cumpos, lg.tick_cumpos)
}

########################################### Look at null and "test" LODs ###################################

#Max LOD score for each comparison and chromosome
max_lods <- aggregate(cbind(S1_v_S2, T1_v_T2,S1_v_T1, S2_v_T2) ~ CHROM, data=lods, FUN=max, na.rm=T)

#What are the genomewide max LOD score for the two null comparisons
max_nulls <- apply(max_lods[,c("S1_v_S2", "T1_v_T2")],2,max)

#Which chromsomes feature LODs > 4 for both of the test comparisons
qtl_chromosomes <- max_lods$CHROM[max_lods$S1_v_T1 > 4 & max_lods$S2_v_T2 >4]

####################################################### Make chromosome plot function ###################################

plot_lg <- function(allele_freqs, #data frame with column for CHROM, POS, and for each pool
                    scaffold_sizes, #data frame where first column is CHROM, second column is chrom length
                    chrom, #which chrom to plot
                    rep, #which rep to plot
                    window_size = 500000, #window size for smooothing
                    step_size = 100000, #step_size for smoothing
                    min_snps = 30, #if a given window has fewer than this many snps, write NA
                    trim_percent = 1, #Show this proportion of SNPs as dots on the plot
                    FUN = function(x) mean(x, na.rm=T), #Function to perform on SNPs in the window
                    plot_arguments = list(), #arguments to pass to plot
                    return_smoothed = F, #If true, return data frame with window smoothed estimates
                    make_plot = T, #whether or not to plot
                    x_axis=F, #Make X axis marks?
                    y_axis=F # Make Y axis marks?
){
  
  require(reshape2)
  FUN <- match.fun(FUN)
  allele_freqs.lg <- allele_freqs[allele_freqs$CHROM==chrom,]
  pos.max <- scaffold_sizes$Length[scaffold_sizes$CHROM == chrom]
  
  #Get window start and stop positions
  start_positions <- seq(0, step_size*ceiling(pos.max/step_size)-window_size, by = step_size)
  if(pos.max %% step_size == 0){
    stop_positions <-  seq(window_size, pos.max , by = step_size)
  }else{
    stop_positions <-  c(seq(window_size, pos.max , by = step_size), pos.max)
  }
  
  #Get window centers based on physical position
  centers.n <- length(start_positions)
  centers <- (stop_positions + start_positions)/2
  smoothed_data <- matrix(NA, nrow=centers.n, ncol=5)
  
  #Now get window smoothed means
  for(i in 1:centers.n){
    window_data <- allele_freqs.lg[allele_freqs.lg$POS > start_positions[i] & allele_freqs.lg$POS <= stop_positions[i], 3:ncol(allele_freqs.lg)]
    pos.window <- allele_freqs.lg$POS[allele_freqs.lg$POS > start_positions[i] & allele_freqs.lg$POS <= stop_positions[i]]
    smoothed_data[i,] <- apply(window_data,2,FUN)
    nsnps <- apply(window_data,2,function(x) sum(!is.na(x)))
    smoothed_data[i,which(nsnps < min_snps)] <- NA
  }
  
  #Put it all together in data frame
  chrom_results <- data.frame("CHROM" = rep(chrom, centers.n),
                              "POS" = centers,
                              smoothed_data)
  colnames(chrom_results)[3:ncol(chrom_results)] <- colnames(allele_freqs.lg)[3:ncol(allele_freqs.lg)]
  
  if(return_smoothed == T){
    return(chrom_results)
  }
  
  #Begin plotting
  if(make_plot == T){
    
    #Get columns that correspond to rep pools
    if(rep == 1){
      pool_columns <- 3:4
    }else if(rep == 2){
      pool_columns <- 5:7
    }
    
    #Get random subset for plotting if trim_percent < 1
    allele_freqs.lg <- allele_freqs.lg[sample(nrow(allele_freqs.lg), size = round(trim_percent*nrow(allele_freqs.lg)), replace = F),]
    
    
    #Turn from wide to long
    allele_freqs.long <- melt(allele_freqs.lg, id.vars=1:2, variable.name = "pool", measure.vars = pool_columns)
    allele_freqs.long <- allele_freqs.long[order(allele_freqs.long$POS),]
    
    #Add colors
    allele_freqs.long$color <- '#D55E00'
    allele_freqs.long$color[grep("t", allele_freqs.long$pool)] <- "#0072B2"
    allele_freqs.long$color[grep("r", allele_freqs.long$pool)] <- '#CC79A2'
    
    #Make plot
    total_plot_arguments <- list(x = allele_freqs.long$POS,
                                 y = allele_freqs.long$value,
                                 col = allele_freqs.long$color,
                                 pch = 3, cex=0.4, ylim=c(0.2,0.8),
                                 xlab = "Position (BP)", ylab = "Frequency Pc-NY21 allele",
                                 xaxt = 'n', yaxt='n')
    total_plot_arguments[names(plot_arguments)] <- plot_arguments
    do.call(plot, total_plot_arguments, quote=T)
    lines(chrom_results$POS, chrom_results[,paste("s", rep, sep="")], col='black', lwd='4.5')
    lines(chrom_results$POS, chrom_results[,paste("s", rep, sep="")], col='#D55E00', lwd='2.5')
    lines(chrom_results$POS, chrom_results[,paste("t", rep, sep="")], col='black', lwd='4.5')
    lines(chrom_results$POS, chrom_results[,paste("t", rep, sep="")], col="#0072B2", lwd='2.5')
    if(rep==2){
      lines(chrom_results$POS, chrom_results[,paste("r", rep, sep="")], col='black', lwd='4.5')
      lines(chrom_results$POS, chrom_results[,paste("r", rep, sep="")], col='#CC79A2', lwd='2.5')
    }
    if(x_axis == T){
      axis_positions <- seq(1,max(allele_freqs.long$POS), by=1000000)
      axis(1, at = axis_positions, labels= 1:length(axis_positions), cex.axis=0.9, padj=0)
    }
    if(y_axis == T){
      axis_positions <- seq(0,1,by=0.10)
      axis(2, at = axis_positions, labels = axis_positions,cex.axis=0.9, las=2, padj=0)
    }
  }
}

#######################################    Look at SNP allele frequency standard deviations in sliding windows     #############################################

freq_sds <- plot_lg(allele_freqs, scaffold_sizes, chrom="Cp4.1LG01", rep=2, step_size=500000,
                    make_plot = F, return_smoothed = T, FUN=function(x) sqrt(var(x, na.rm=T)))
for(chrom in unique(allele_freqs$CHROM)){
  freq_sds <- rbind(freq_sds, 
                    plot_lg(allele_freqs, scaffold_sizes, chrom=chrom, rep=2, step_size=500000,
                            make_plot = F, return_smoothed = T, FUN=function(x) sqrt(var(x, na.rm=T))))
}
apply(freq_sds[3:ncol(freq_sds)],2,summary)

apply(allele_freqs[,3:ncol(allele_freqs)],2,function(x) sqrt(var(x,na.rm=T)))

#######################################    Make plots for each chromosome and rep     #############################################

#First make plots
for(rep in 1:2){
  
  pdf(paste("plots/allele_freqs_rep",rep, ".pdf", sep=''), height=10, width=8)
  
  old.par <- par(no.readonly = T)
  m <- rbind(c(1,2,2,3,3,4,4,5,5),
             c(1,2,2,3,3,4,4,5,5),
             c(1,6,6,7,7,8,8,9,9),
             c(1,6,6,7,7,8,8,9,9),
             c(1,10,10,11,11,12,12,13,13),
             c(1,10,10,11,11,12,12,13,13),
             c(1,14,14,15,15,16,16,17,17),
             c(1,14,14,15,15,16,16,17,17),
             c(1,18,18,19,19,20,20,21,21),
             c(1,18,18,19,19,20,20,21,21),
             c(1,22,22,22,22,22,22,22,22))
  layout(m)
  
  par(oma=c(1,1,1,1), mar=c(1,1,1,1))
  
  #Empty plot where some labels will go
  plot(0,xlim=c(0,10),ylim=c(0,10),type='n', xaxt='n',yaxt='n',bty='n', xlab='', ylab='')
  text(3,5.5, "Pc-NY21 Allele Frequency", srt=90, cex=1.5)
  
  #Loop thru chromosomes
  par(mar=c(2,0.5,2,0.5))
  for(chrom in as.character(unique(allele_freqs$CHROM))){
  
    chrom_number <- as.integer(gsub("Cp4.1LG", "", chrom))
    
    #Plot x and y axes?
    y_axis <- F
    if(chrom %in% c("Cp4.1LG01", "Cp4.1LG05", "Cp4.1LG09", "Cp4.1LG13", "Cp4.1LG17")){
      y_axis <- T
    }
    plot_lg(allele_freqs, scaffold_sizes, chrom, rep, x_axis=T, y_axis=y_axis, trim_percent = 0.25)
    
    #Add chromosome label
    plot_center <- grconvertX(0.5, "npc", "user")
    left_bound <- grconvertX(0.3, "npc", "user")
    right_bound <- grconvertX(0.7, "npc", "user")
    plot_top <- grconvertY(1, "npc", "user")
    rect(left_bound, 0.72, right_bound, plot_top, col="white", border="black")
    text(plot_center, 0.77, paste("CHROM", chrom_number))
  }
  
  #Legend at bottom
  par(mar=c(0,0.5,2,0.5))
  plot(0, xlim=c(0,10), ylim=c(0,10),type='n', xaxt='n', yaxt='n', bty='n', ylab = 'n', xlab = 'n')
  text(5,11, "Position (Mb)", cex=1.5, xpd=NA)
  if(rep == 1){
    legend("bottom",
           fill = c('#D55E00', "#0072B2"),
           legend = c("SUS Pool", "RES pool"),
           bty = 'n', ncol=3, cex=1.75)
  }else if(rep == 2){
    legend("bottom",
           fill = c('#D55E00', "#0072B2", '#CC79A2'),
           legend = c("SUS Pool", "RES pool", "RAN pool"),
           bty = 'n', ncol=3, cex=2, x.intersp=1.5)
    
  }
  par(old.par)
  dev.off()
}

############################################    Make table with QTL lods and max AF diff     #############################################

#Get the location and magnitude of the maximum smoothed allele frequency difference per chromosome
chroms <- unique(allele_freqs$CHROM)
max_diff <- data.frame("CHROM" = chroms,
                       "Pos_R1" = NA,
                       "DeltaAF_R1" = NA,
                       "Pos_R2" = NA,
                       "DeltaAF_R2" = NA)
for(chrom in chroms){
  for(rep in 1:2){
    af_smooth <- plot_lg(allele_freqs, scaffold_sizes, chrom, rep, make_plot = F, return_smoothed = T)
    diff <- af_smooth[,paste("t",rep,sep="")] - af_smooth[,paste("s",rep,sep="")]
    diff_coord <- which(abs(diff)==max(abs(diff), na.rm = T))
    pull <- c(af_smooth$POS[diff_coord], diff[diff_coord])
    max_diff[max_diff$CHROM == chrom, c(rep*2, rep*2+1)] <- pull
  }
}

write.csv(max_diff, "tables/max_smoothed_deltaAF.csv", quote=F, row.names = F)

#Make final table, combining multipool output with allele frequency differences for QTL chromosomes

#Function for converting to Mb and rounding to 2 sig digs
Bp_to_Mb <- function(x){
  return(round(x/1000000,2))
}

#Function for listing Rep 1 and Rep 2 values together like:
# "XXXX / YYYYY"
combine_reps <- function(rep1, rep2){
  return(paste(rep1, "/", rep2, sep=" "))
}

#Function for turning peak and CI start and end into something like:
# "Pos (Start - End)"
addCI <- function(peak, start, end){
  return(paste(peak, "(", start, "-", end, ")", sep=""))
}
chrom_summary <- data.frame("Chrom" = max_diff$CHROM,
                            "MaxLOD" = combine_reps(max_lods$S1_v_T1, max_lods$S2_v_T2),
                            "Pos_CI" = combine_reps(addCI(multi_peaks$Rep1_peak, multi_peaks$Rep1_start, multi_peaks$Rep1_end),
                                                    addCI(multi_peaks$Rep2_peak, multi_peaks$Rep2_start, multi_peaks$Rep2_end)),
                            "MaxDeltaAF" = combine_reps(max_diff$DeltaAF_R1, max_diff$DeltaAF_R2),
                            "Pos" = combine_reps(max_diff$Pos_R1, max_diff$Pos_R2))
chrom_summary.round <- data.frame("Chrom" = max_diff$CHROM,
                          "MaxLOD" = combine_reps(max_lods$S1_v_T1, max_lods$S2_v_T2),
                          "Pos_CI" = combine_reps(addCI(Bp_to_Mb(multi_peaks$Rep1_peak), Bp_to_Mb(multi_peaks$Rep1_start), Bp_to_Mb(multi_peaks$Rep1_end)),
                                                  addCI(Bp_to_Mb(multi_peaks$Rep2_peak), Bp_to_Mb(multi_peaks$Rep2_start), Bp_to_Mb(multi_peaks$Rep2_end))),
                          "MaxDeltaAF" = combine_reps(round(max_diff$DeltaAF_R1,2), round(max_diff$DeltaAF_R2,2)),
                          "Pos" = combine_reps(Bp_to_Mb(max_diff$Pos_R1), Bp_to_Mb(max_diff$Pos_R2)))
qtl_summary <- chrom_summary[chrom_summary$Chrom %in% qtl_chromosomes,]
qtl_summary.round <- chrom_summary.round[chrom_summary.round$Chrom %in% qtl_chromosomes,]


write.csv(qtl_summary, "tables/BSA_results.csv", quote=F, row.names = F)
write.csv(qtl_summary.round, "tables/BSA_results_round.csv", quote=F, row.names = F)

################### Make multipane plot with LOD scores and allele freqs on chroms 4,5,8,16 ######################################

pdf("plots/BSA_results.pdf", width=6.85089, height=5.872191)
old.par <- par(no.readonly = T)
par(oma=c(0.5,0,2,1))

m <- rbind(c(1,2,3,4,5,6,12),
           c(1,7,8,9,10,11,12 ),
           13,
           14,
           14)
layout(m)

par(mar=c(0.25,0,0.25,0))
#Empty plot where some labels will go
plot(0,xlim=c(0,10),ylim=c(0,10),type='n', xaxt='n',yaxt='n',bty='n', xlab='', ylab='')
text(5,5, "Pc-NY21 Allele Frequency", srt=90, cex=1.2)

#Add letters
par(xpd=NA)
text(2.5,10.5, "A", cex=2)
text(2.5,-5.5, "B", cex=2)
par(xpd=F)

par(mar=c(0.25,1,0.25,1))
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG04", 1, y_axis=T, trim_percent = 0.25)
mtext('Chrom 4',side=3, line=0.25, cex=0.75)
mtext('Rep 1', side=2, line=2.5, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG05", 1, trim_percent = 0.25)
mtext('Chrom 5',side=3, line=0.25, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG08", 1, trim_percent = 0.25)
mtext('Chrom 8',side=3, line=0.25, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG12", 1, trim_percent = 0.25)
mtext('Chrom 12',side=3, line=0.25, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG16", 1, trim_percent = 0.25)
mtext('Chrom 16',side=3, line=0.25, cex=0.75)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG04", 2, x_axis=T, y_axis=T, trim_percent=0.25)
mtext('Rep 2', side=2, line=2.5, cex=0.75)
mtext('Position (Mb)', side=1, line=2.5, cex=0.75)
y_end <- grconvertY(0, "npc", "ndc")
chr4_end1 <- grconvertX(0, "npc", "ndc")
chr4_end2 <- grconvertX(1, "npc", "ndc")

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG05", 2, x_axis=T, trim_percent=0.25)
chr5_end1 <- grconvertX(0, "npc", "ndc")
chr5_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2.5, cex=0.75)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG08", 2, x_axis=T, trim_percent=0.25)
chr8_end1 <- grconvertX(0, "npc", "ndc")
chr8_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2.5, cex=0.75)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG12", 2, x_axis=T, trim_percent=0.25)
chr12_end1 <- grconvertX(0, "npc", "ndc")
chr12_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2.5, cex=0.75)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG16", 2, x_axis=T, trim_percent=0.25)
chr16_end1 <- grconvertX(0, "npc", "ndc")
chr16_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2.5, cex=0.75)

#Legend plot
par(mar=c(0.25,0,0.25,0))
plot(0,xlim=c(0,10),ylim=c(0,10),type='n', xaxt='n',yaxt='n',bty='n', xlab='', ylab='')
legend(0,7,
       fill = c('#D55E00', "#0072B2", '#CC79A2'), 
       legend = c("SUS Pool", "RES pool", "RAN pool"), 
       bty = 'n',cex=1.25)

#Add blank plot in between
plot(0,xlim=c(0,10),ylim=c(0,10),type='n', xaxt='n',yaxt='n',bty='n', xlab='', ylab='')

#Add LOD scores Plot
par(mar=c(4,4,1,0.25))
plot(x=0, type="n", xaxt="n",
     xlim=range(cumpos),
     ylim = range(c(lods$S1_v_T1, lods$S2_v_T2), na.rm=T),
     xlab = "", ylab="",
     cex.axis=1, las=2)
mtext("LOD", side=2, line=1.75, cex=0.75)
mtext("Genomic position", side=1, line=2, cex=0.75)

for(lg in lgs.unique){
  lines(x=cumpos[lgs==lg],
        y=lods$S1_v_T1[lgs==lg],
        col="tan1", lwd=1.25)
  lines(x=cumpos[lgs==lg],
        y=lods$S2_v_T2[lgs==lg],
        col="lightskyblue4", lwd=1.25)
}

for(i in 1:length(lgs.cumsum)){
  abline(v=lgs.cumsum[i], col='black', lty=2,lwd=1.25)
}

legend('topright', fill = c("tan1", "lightskyblue4"), legend = c("Rep 1", "Rep 2"), bg="white")

axis(side = 1, at = lgs.centers[,2],
     labels = 1:20,
     tick = F)

y_start <- grconvertY(1,"npc","ndc")
chr4_start1 <- grconvertX(lgs.cumsum[4], "user", "ndc")
chr4_start2 <- grconvertX(lgs.cumsum[5], "user", "ndc")
chr5_start1 <- grconvertX(lgs.cumsum[5], "user", "ndc")
chr5_start2 <- grconvertX(lgs.cumsum[6], "user", "ndc")
chr8_start1 <- grconvertX(lgs.cumsum[8], "user", "ndc")
chr8_start2 <- grconvertX(lgs.cumsum[9], "user", "ndc")
chr12_start1 <- grconvertX(lgs.cumsum[12], "user", "ndc")
chr12_start2 <- grconvertX(lgs.cumsum[13], "user", "ndc")
chr16_start1 <- grconvertX(lgs.cumsum[16], "user", "ndc")
chr16_start2 <- grconvertX(lgs.cumsum[17], "user", "ndc")
abline(h=4, col='red')

pushViewport(viewport())
grid.lines(x=c(chr4_start1,chr4_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr4_start2,chr4_end2), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr5_start1,chr5_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr5_start2,chr5_end2), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr8_start1,chr8_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr8_start2,chr8_end2), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr12_start1,chr12_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr12_start2,chr12_end2), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr16_start1,chr16_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr16_start2,chr16_end2), y=c(y_start,y_end), gp=gpar(lty=2))

axis(side=1, at =tick_cumpos, line=0,
     labels = FALSE, tick = T, cex.axis=0.5, lwd.ticks = 0.75)

dev.off()

