setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/BSA/")
library(reshape2)

#########################################################   Import data and load scripts  #############################################################
source("../scripts/Pcap-QTL-Mapping/BSA_analysis/Random_functions_for_sliding_window_analyses.R")

counts <- read.table("data/counts_filtered.txt", stringsAsFactors = F, header=T)
lods <- read.table("data/multipool_lods.txt", header=T)
scaffold_sizes <- read.table("data/scaffold_sizes.txt")

####################################################### Process data prior to making plot ###################################

#Add column names to scaffold sizes and remove LG 00
colnames(scaffold_sizes) <- c("CHROM", "Length")
scaffold_sizes <- scaffold_sizes[scaffold_sizes$CHROM!="Cp4.1LG00",]
scaffold_sizes$CHROM <- droplevels(scaffold_sizes$CHROM)

#Turn allele counts into freqs
allele_freqs <- collapse_columns(counts = counts, data_columns = 3:ncol(counts), stat_type = "frequency",
                                 rd_filter = c(0,1), rd_filter_type = "quantile", freq_filter = 1)

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
####################################################### Make chromosome plot function ###################################

plot_lg <- function(allele_freqs, #data frame with column for CHROM, POS, and for each pool
                    scaffold_sizes, #data frame where first column is CHROM, second column is chrom length
                    chrom, #which chrom to plot
                    rep, #which rep to plot
                    window_size = 500000, #window size for smooothing
                    step_size = 100000, #step_size for smoothing
                    min_snps = 30, #if a given window has fewer than this many snps, write NA
                    trim_percent = 1, #Show this proportion of SNPs as dots on the plot
                    plot_arguments = list(), #arguments to pass to plot
                    return_smoothed = F, #If true, return data frame with window smoothed estimates
                    make_plot = T, #whether or not to plot
                    x_axis=F, #Make X axis marks?
                    y_axis=F # Make Y axis marks?
){
  
  require(reshape2)
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
    smoothed_data[i,] <- apply(window_data,2,mean)
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
    lines(chrom_results$POS, chrom_results[,paste("s", rep, sep="")], col='black', lwd='6')
    lines(chrom_results$POS, chrom_results[,paste("s", rep, sep="")], col='#D55E00', lwd='4')
    lines(chrom_results$POS, chrom_results[,paste("t", rep, sep="")], col='black', lwd='6')
    lines(chrom_results$POS, chrom_results[,paste("t", rep, sep="")], col="#0072B2", lwd='4')
    if(rep==2){
      lines(chrom_results$POS, chrom_results[,paste("r", rep, sep="")], col='black', lwd='6')
      lines(chrom_results$POS, chrom_results[,paste("r", rep, sep="")], col='#CC79A2', lwd='4')
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

################################## Make plots for each chromosome and rep and pull max AF difference ######################################

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
    plot_lg(allele_freqs, scaffold_sizes, chrom, rep, x_axis=T, y_axis=y_axis)
    
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

################### Make multipane plot with LOD scores and allele freqs on chroms 4,5,8,16 ######################################

pdf("plots/BSA_results.pdf", width=7, height=6)
old.par <- par(no.readonly = T)
par(oma=c(0.5,0,2,1))

m <- rbind(c(1,2,3,4,5,10),
           c(1,6,7,8,9,10 ),
           c(11,11,11,11,11,11),
           c(12,12,12,12,12,12),
           c(12,12,12,12,12,12))
layout(m)

par(mar=c(0.25,0,0.25,0))
#Empty plot where some labels will go
plot(0,xlim=c(0,10),ylim=c(0,10),type='n', xaxt='n',yaxt='n',bty='n', xlab='', ylab='')
text(6,5, "Pc-NY21 Allele Frequency", srt=90, cex=1.2)

par(mar=c(0.25,1,0.25,1))
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG04", 1, y_axis=T, trim_percent = 0.25)
mtext('Chrom 4',side=3, line=0.25, cex=0.75)
mtext('Rep 1', side=2, line=2.5, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG05", 1, trim_percent = 0.25)
mtext('Chrom 5',side=3, line=0.25, cex=0.75)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG08", 1, trim_percent = 0.25)
mtext('Chrom 8',side=3, line=0.25, cex=0.75)
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
     ylim = range(c(lods$S1_v_T1, lods$S2_v_T2)),
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
grid.lines(x=c(chr16_start1,chr16_end1), y=c(y_start,y_end), gp=gpar(lty=2))
grid.lines(x=c(chr16_start2,chr16_end2), y=c(y_start,y_end), gp=gpar(lty=2))

axis(side=1, at =tick_cumpos, line=0,
     labels = FALSE, tick = T, cex.axis=0.5, lwd.ticks = 0.75)

dev.off()

