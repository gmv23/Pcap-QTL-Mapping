setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/synteny/")
library(circlize)

#Read in input files: file with info on matching syntenic blocks
#And files with length of chromosomes for each of Cpepo and Cmoschata
#And file with QTL coordinates
syn <- read.table("squash.collinearity.flat", header=T)
Cm_sizes <- read.table("Cm_chrom_lengths.txt")
Cp_sizes <- read.table("Cp_chrom_lengths.txt")
qtl <- read.csv("qtl.csv", stringsAsFactors = F)

#Clean up chrom sizes tables and rename chromosomes to integers and add species Id
colnames(Cm_sizes) <- c("Chrom", "Size")
Cm_sizes$Chrom <- as.integer(gsub("Cmo_Chr", "", Cm_sizes$Chrom))
Cm_sizes$ChromID <- paste("Cm", Cm_sizes$Chrom, sep="")
colnames(Cp_sizes) <- c("Chrom", "Size")
Cp_sizes$Chrom <- as.integer(gsub("Cp4.1LG", "", Cp_sizes$Chrom))
Cp_sizes$ChromID <- paste("Cp", Cp_sizes$Chrom, sep="")

#Subset syn blocks by statistical significance to only show most important?
syn <- syn[syn$Score > 400,]


pdf("plots/synteny_plot.pdf")
m <- rbind(c(1,1,2,2),c(4,3,3,4))
layout(m)

for(Cm_chrom in c(4,11,14)){
  plot_chroms <- rbind(Cm_sizes[Cm_sizes$Chrom==Cm_chrom,], Cp_sizes[Cp_sizes$Chrom!=0,])
  #Add y values with range (0,1)
  plot_chroms$y <- 1
  #Add 0 for each chromosome coordinate so that we have effectively a range for each factor
  plot_chroms0 <- plot_chroms
  plot_chroms0$Size <- 0
  plot_chroms0$y <- 0
  plot_chroms <- rbind(plot_chroms, plot_chroms0)
  #Add background color
  plot_chroms$bg <- "gray"
  plot_chroms$bg[grep("Cm",plot_chroms$ChromID)] <- "blue"
  
  #Reorder factors
  plot_chroms$ChromID <- as.factor(plot_chroms$ChromID)
  factor_order <- c(paste("Cm",Cm_chrom, sep=""), paste("Cp",1:20, sep=""))
  plot_chroms$ChromID <- factor(plot_chroms$ChromID, levels=factor_order)
  
  #Initalize circle
  circos.par("track.height" = 0.1, "start.degree" = 100)
  circos.initialize(factors=plot_chroms$ChromID, x=plot_chroms$Size)
  
  #Add track with labels
  circos.track(factors = plot_chroms$ChromID, y = plot_chroms$y, bg.col=plot_chroms$bg,
               panel.fun = function(x,y){
                 circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), CELL_META$sector.index)
               })
  
  #Draw QTL regions
  qtl_plot <- qtl[qtl$Chrom %in% plot_chroms$ChromID,]
  for(i in 1:nrow(qtl_plot)){
    circos.rect(qtl_plot$Start[i], -0.2, qtl_plot$Stop[i], 1.2, sector.index=qtl_plot$Chrom[i], col='red', border=NA)
  }
  
  #Draw links
  syn_plot <- syn[syn$Chrom_Cm == Cm_chrom & syn$Chrom_Cp != 0,]
  syn_plot$Chrom_Cm <- paste("Cm", syn_plot$Chrom_Cm, sep="")
  syn_plot$Chrom_Cp <- paste("Cp", syn_plot$Chrom_Cp, sep="")
  
  #Pull out regions of Cm QTl we are comparing to color syntenic blocks in that region red
  qtl_min <- qtl_plot[qtl_plot$Chrom==paste("Cm",Cm_chrom,sep=""), "Start"]
  qtl_max <- qtl_plot[qtl_plot$Chrom==paste("Cm",Cm_chrom,sep=""), "Stop"]
  
  syn_plot_qtl <- syn_plot[((syn_plot$Start_Cm > qtl_min & syn_plot$Start_Cm < qtl_max) | 
             (syn_plot$Stop_Cm > qtl_min & syn_plot$Stop_Cm < qtl_max)),] 
  syn_plot <- syn_plot[-which(((syn_plot$Start_Cm > qtl_min & syn_plot$Start_Cm < qtl_max) | 
                          (syn_plot$Stop_Cm > qtl_min & syn_plot$Stop_Cm < qtl_max))),] 
  
  for(i in 1:nrow(syn_plot)){
    circos.link(sector.index1=syn_plot$Chrom_Cm[i], c(syn_plot$Start_Cm[i], syn_plot$Stop_Cm[i]),
                sector.index2=syn_plot$Chrom_Cp[i], c(syn_plot$Start_Cp[i], syn_plot$Stop_Cp[i]),
                col=adjustcolor('gray', alpha.f=0.6))
  }
  for(i in 1:nrow(syn_plot_qtl)){
    circos.link(sector.index1=syn_plot_qtl$Chrom_Cm[i], c(syn_plot_qtl$Start_Cm[i], syn_plot_qtl$Stop_Cm[i]),
                sector.index2=syn_plot_qtl$Chrom_Cp[i], c(syn_plot_qtl$Start_Cp[i], syn_plot_qtl$Stop_Cp[i]),
                col=adjustcolor('red', alpha.f=0.6))
  }
  circos.clear()
}
dev.off()


