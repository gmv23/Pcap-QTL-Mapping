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
                    window_size = 1000000, #window size for smooothing
                    step_size = 100000, #step_size for smoothing
                    min_snps = 30, #if a given window has fewer than this many snps, write NA
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
                                 pch = 3, cex=0.6, ylim=c(0.2,0.8),
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
      axis(1, at = axis_positions, labels= 1:length(axis_positions), cex.axis=0.6)
    }
    if(y_axis == T){
      axis_positions <- seq(0,1,by=0.10)
      axis(2, at = axis_positions, labels = axis_positions,cex.axis=0.6, las=2)
    }
  }
}

#Make a plot for each chromosome
for(chrom in as.character(unique(allele_freqs$CHROM))){
  for(rep in 1:2){
    pdf(paste("plots/", "Rep", rep, "_", chrom, ".pdf", sep=''))
    plot_lg(allele_freqs, scaffold_sizes, chrom, rep)
    dev.off()
  }
}

#################################################### Make multipane plot ######################################

pdf("BSA_results.pdf", width=7, height=5)
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
text(6,5, "Pc-NY21 Allele Frequency", srt=90)

par(mar=c(0.25,0.75,0.75,1))
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG04", 1, y_axis=T)
mtext('Chrom 4',side=3, line=0.25, cex=0.6)
mtext('Rep 1', side=2, line=2, cex=0.6)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG05", 1)
mtext('Chrom 5',side=3, line=0.25, cex=0.6)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG08", 1)
mtext('Chrom 8',side=3, line=0.25, cex=0.6)
plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG16", 1)
mtext('Chrom 16',side=3, line=0.25, cex=0.6)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG04", 2, x_axis=T, y_axis=T)
mtext('Rep 2', side=2, line=2, cex=0.6)
mtext('Position (Mb)', side=1, line=2, cex=0.6)
y_end <- grconvertY(0, "npc", "ndc")
chr4_end1 <- grconvertX(0, "npc", "ndc")
chr4_end2 <- grconvertX(1, "npc", "ndc")

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG05", 2, x_axis=T)
chr5_end1 <- grconvertX(0, "npc", "ndc")
chr5_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2, cex=0.6)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG08", 2, x_axis=T)
chr8_end1 <- grconvertX(0, "npc", "ndc")
chr8_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2, cex=0.6)

plot_lg(allele_freqs, scaffold_sizes, "Cp4.1LG16", 2, x_axis=T)
chr16_end1 <- grconvertX(0, "npc", "ndc")
chr16_end2 <- grconvertX(1, "npc", "ndc")
mtext('Position (Mb)', side=1, line=2, cex=0.6)

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
     xlab = "Genomic position",
     cex.axis=0.8, las=2)
mtext("LOD", side=2, line=2)

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




































############### Look how plot altogether seems

pdf("plots/BSA_results.pdf")
old.par <- par(no.readonly = T)
par(mar=c(1,1,1,1), oma=c(4,6,1,1))
m <- rbind(c(1,3,5,7),
           c(2,4,6,8),
           c(9,9,9,9))
layout(m)
for(lg in lgs.unique[c(4,5,8,16)]){ #loop thru lgs
  allele_freqs.lg <- allele_freqs[allele_freqs$CHROM==lg,]
  
  pos.max <- scaffold_sizes[scaffold_sizes[,1] == lg,2]
  
  #Window smoothed estimates
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
  
  for(i in 1:centers.n){
    window_data <- allele_freqs.lg[allele_freqs.lg$POS > start_positions[i] & allele_freqs.lg$POS <= stop_positions[i], 3:ncol(allele_freqs.lg)]
    pos.window <- allele_freqs.lg$POS[allele_freqs.lg$POS > start_positions[i] & allele_freqs.lg$POS <= stop_positions[i]]
    smoothed_data[i,] <- apply(window_data,2,mean)
    nsnps <- apply(window_data,2,function(x) sum(!is.na(x)))
    smoothed_data[i,which(nsnps < min_snps)] <- NA
  }
  
  chrom_results <- data.frame("CHROM" = rep(lg, centers.n),
                              "POS" = centers,
                              smoothed_data)
  colnames(chrom_results)[3:ncol(chrom_results)] <- colnames(allele_freqs.lg)[3:ncol(allele_freqs.lg)]
  
  #Turn from wide to long
  rep1_freqs <- melt(allele_freqs.lg[,1:4], id.vars=1:2, variable.name = "pool")
  rep1_freqs <- rep1_freqs[order(rep1_freqs$POS),]
  rep2_freqs <- melt(allele_freqs.lg[,c(1:2,5:7)], id.vars=1:2, variable.name = "pool")
  rep2_freqs <- rep2_freqs[order(rep2_freqs$POS),]
  
  #Assign color to pools
  rep1_freqs$color <- '#D55E00'
  rep1_freqs$color[rep1_freqs$pool == "t1"] <- "#0072B2"
  rep2_freqs$color <- '#D55E00'
  rep2_freqs$color[rep2_freqs$pool == "t2"] <- "#0072B2"
  rep2_freqs$color[rep2_freqs$pool == "r2"] <- '#CC79A2'
  
  #Now start making plots
  yaxt='n'
  if(lg=="Cp4.1LG04"){
    yaxt='s'
  }
  plot(rep1_freqs$POS, rep1_freqs$value, pch=3, cex=0.6,col=rep1_freqs$color,
       ylim=c(0.2,0.8), xaxt='n', yaxt=yaxt, main=lg, ylab='Pc-NY21 Allele Frequency \n Rep 1')
  lines(chrom_results$POS, chrom_results$s1, col='black', lwd='6')
  lines(chrom_results$POS, chrom_results$s1, col='#D55E00', lwd='4')
  lines(chrom_results$POS, chrom_results$t1, col='black', lwd='6')
  lines(chrom_results$POS, chrom_results$t1, col="#0072B2", lwd='4')

  plot(rep2_freqs$POS, rep2_freqs$value, pch=3, cex=0.5,col=rep2_freqs$color, 
       ylim=c(0.2,0.8), yaxt=yaxt, ylab='Pc-NY21 Allele Frequency \n Rep 2')
  lines(chrom_results$POS, chrom_results$s2, col='black', lwd='6')
  lines(chrom_results$POS, chrom_results$s2, col='#D55E00', lwd='4')
  lines(chrom_results$POS, chrom_results$t2, col='black', lwd='6')
  lines(chrom_results$POS, chrom_results$t2, col="#0072B2", lwd='4')
  lines(chrom_results$POS, chrom_results$r2, col='black', lwd='6')
  lines(chrom_results$POS, chrom_results$r2, col='#CC79A2', lwd='4')
}


par(mar=c(1,1,3,1),xpd='T')
plot(x=0, type="n", xaxt="n",
     xlim=range(cumpos),
     ylim = range(c(lods$S1_v_T1, lods$S2_v_T2)),
     xlab = "Genomic position", ylab="LOD",
     cex.lab=1.3,
     cex.main = 1.5)

#add lines
for(lg in lgs.unique){
  lines(x=cumpos[lgs==lg],
        y=lods$S1_v_T1[lgs==lg],
        col='orange')
  lines(x=cumpos[lgs==lg],
        y=lods$S2_v_T2[lgs==lg],
        col='blue')
}
for(i in 1:length(lgs.cumsum)){
  abline(v=lgs.cumsum[i], col='gray', lty=2,xpd=F)
}

axis(side = 1, at = lgs.centers[,2],
         labels = 1:20,
         tick = F, line=0.5, cex.axis=1)

abline(h=4,xpd=F)

par(old.par)
dev.off()

