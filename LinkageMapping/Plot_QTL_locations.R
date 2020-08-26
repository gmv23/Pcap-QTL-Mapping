#Plot CIs and peak position of QTLs as detected by BSA Rep 1, BSA Rep 2, and linkage mapping


setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/LinkageMapping/")

scaffold_sizes <- read.table("../BSA/data/scaffold_sizes.txt")
multi_peaks <- read.csv("../BSA/data/CIs_and_peaks.txt")
snps <- read.csv("../genetic_map/tables/snps_filter_prune.csv")
lm_peaks <- read.csv("tables/qtl_summary_lm.csv", stringsAsFactors = F)
gbs_markers <- read.csv("tables/GBS_markers_print.csv")

colnames(scaffold_sizes) <- c("CHROM", "Length")
scaffold_sizes <- scaffold_sizes[scaffold_sizes$CHROM!="Cp4.1LG00",]
scaffold_sizes$CHROM <- droplevels(scaffold_sizes$CHROM)

make_lg_integer <- function(x){
  return(as.integer(gsub("Cp4.1LG","",x)))
}

#Make chromosome names integers
scaffold_sizes$CHROM <- sapply(scaffold_sizes$CHROM, make_lg_integer)
snps$CHROM <- sapply(snps$CHROM, make_lg_integer)
multi_peaks$CHROM <- sapply(multi_peaks$CHROM, make_lg_integer)

#Subset just QTL chromosomes in multi_peaks and divide into two data frames for each rep
multi_peaks <- multi_peaks[multi_peaks$CHROM %in% c(4,5,8,12,16),]

#Rename columns

#Subset scaffold_sizes for just chromosomes that have a QTL
chroms.qtl <- c(4,5,8,12,16,19)
scaffold_sizes <- scaffold_sizes[scaffold_sizes$CHROM %in% chroms.qtl,]

y_coords <- seq(100,0,length.out=6)


#pdf("plots/QTL_positions.pdf", height=5, width=7)
jpeg("plots/QTL_positions.jpeg", height = 5, width=7, units="in", res=500)
old.par <- par(no.readonly = T)
par(mar=c(5,6,0,7), xpd=NA)
plot(0, type='n',xlim=c(0,max(scaffold_sizes)),ylim=c(0,110),
     bty='n', xaxt='n', yaxt='n',ylab='',xlab='Position (Mb)')

for(i in 1:length(chroms.qtl)){
  chrom <- chroms.qtl[i]
  
  #Draw line for chromosome
  y_coord <- y_coords[i]
  segments(x0=0,x1=scaffold_sizes$Length[scaffold_sizes$CHROM==chrom],
         y0=y_coord, y1=y_coord)
  
  #Add marker tick marks
  markers <- snps$BP[snps$CHROM==chrom]
  for(marker in markers){
    segments(x0=marker,x1=marker,
           y0=y_coord-1,y1=y_coord+1)
  }
  
  #Add GBS marker tagging BSA QTL
  if(chrom %in% gbs_markers$CHROM){
    points(gbs_markers$GBS_BP[gbs_markers$CHROM==chrom],
           y_coord,
           pch=8,cex=1)
  }
  
  #Begin drawing confidence intervals
  interval_number <- 0
  
  if(chrom %in% multi_peaks$CHROM){
    #Draw Rep 1 BSA interval
    interval_number <- interval_number + 1
    interval_pos <- y_coord + (interval_number * 3)
    segments(x0 = multi_peaks$Rep1_start[multi_peaks$CHROM==chrom],
             x1 = multi_peaks$Rep1_end[multi_peaks$CHROM==chrom],
             y0 = interval_pos, y1= interval_pos,
             lwd=3, col="orange")
    points(multi_peaks$Rep1_peak[multi_peaks$CHROM==chrom], interval_pos,
        pch=21, col='black', bg='white', cex=1.5)
    #Draw Rep 2 BSA interval
    interval_number <- interval_number + 1
    interval_pos <- y_coord + (interval_number * 3)
    segments(x0 = multi_peaks$Rep2_start[multi_peaks$CHROM==chrom],
             x1 = multi_peaks$Rep2_end[multi_peaks$CHROM==chrom],
             y0 = interval_pos, y1= interval_pos,
             lwd=3, col="red")
    points(multi_peaks$Rep2_peak[multi_peaks$CHROM==chrom], interval_pos,
           pch=21, col='black', bg='white', cex=1.5)
  }
  
  if(chrom %in% lm_peaks$Chrom){
    interval_number <- interval_number + 1
    interval_pos <- y_coord + (interval_number*3)
    
    cis <- lm_peaks$CI[lm_peaks$Chrom==chrom]
    #Some QTL have disjoint CIs so we have to loop through them
    cis <- unlist(strsplit(cis,";"))
    for(ci in cis){
      bounds <- as.integer(unlist(strsplit(ci,"-")))
      segments(x0 = bounds[1],
               x1 = bounds[2],
               y0 = interval_pos, y1= interval_pos,
               lwd=3, col="purple")
    }
    points(lm_peaks$Peak[lm_peaks$Chrom==chrom], interval_pos,
           pch=21, col='black', bg='white', cex=1.5)
  }
  
  #Add x lab
  axis(side=1,at=seq(0, max(scaffold_sizes$Length), by=1000000),
       labels=seq(0, max(scaffold_sizes$Length), by=1000000)/1000000)
  
  #Add chromosome labels
  text(rep(-2000000,6), y_coords,
       labels=paste("Chromosome\n",chroms.qtl))
  
  #Add legend
  legend(max(scaffold_sizes$Length) - 1200000, 80, bty='n',
         legend = c("BSA-Seq Rep1", "BSA-Seq Rep2", "MQM"),
         col = c("orange", "red", "purple"), lwd=3)
  
  legend(max(scaffold_sizes$Length) - 1200000, 50, bty='n',
         legend = "BSA-Seq QTL\ntagging marker",
         pch=8)
  
}
par(old.par)
dev.off()
