setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/genetic_map/")
library(qtl)
library(ASMap)

#Import data
squash <- read.cross(file="tables/f2_geno_rqtl.csv", format='csv')
snps <- read.csv("tables/snps_filter_prune.txt")

########################      LOOK AT RECOMBINATION FRACTIONS BETWEEN MARKERS    ###################

#Look at map without any genetic distance info (just # of markers per chromosome)
plotMap(squash)

#Look at recombination fractions and LD scores between all marker pairs
squash <- est.rf(squash)

#Plot recombination fractors
pdf("plots/RF.pdf")
plotRF(squash)
dev.off()

#Plot RFs for each chromosome
for(i in 1:20){
  pdf(paste('plots/RF_chrom', i, '.pdf', sep=''))
  plotRF(squash, chr=i)
  dev.off()
}

#Make map and plot
newmap <- est.map(squash, error.prob=0.003, map.function="kosambi")
plotMap(newmap)
#########################          FIX MARKER PLACEMENT PROBLEMS        ################################

#There appears to be quite a few errors in marker placement order:
# a marker on chrom 1 -> looks like it belongs on chrom 12
# a block of 2 markers on chrom 2 -> looks like they belong to group of markers on chroms 9 and 10
# Something weird on chrom 4 --- maybe big inversion
# A marker on chrom 9 that looks like it goes on the beginning of chrom 15
# Two markers on chrom 10 that belong on chrom 9 along with the 2 markers on chrom 2
# Another marker on chrom 10 that belongs on chrom 13
# Maybe some sort of big inversion on chrom 17

#Move markers to correct LG; dont worry about placement as all marker orders will be re-estimated
#Chrom 1
plotMap(newmap, chr = 1, show.marker.names = T)
plotRF(squash, chr = c(1,12))
squash <- movemarker(squash, "1_5829601", 12)
plotRF(squash, chr = c(1,12))

#Chrom 9
plotMap(newmap, chr = 9, show.marker.names = T)
plotRF(squash, chr = c(9,15))
squash <- movemarker(squash, "9_5600527", 15) 
plotRF(squash, chr=c(9,15))

#Chrom 10 -- markers on chromosome 9:9479477, 9557605
plotRF(squash, chr = c(10,9))
plotMap(newmap, chr = 10, show.marker.names = T)
squash <- movemarker(squash, "10_9479477", 9) 
squash <- movemarker(squash, "10_9557605", 9) 

plotRF(squash, chr = c(9,10,13))

#Chrom 2 -- markers on chromosome 9: 2_8175080, 8224962
plotRF(squash, chr = c(2,9))
plotMap(newmap, chr = 2, show.marker.names = T)
squash <- movemarker(squash, "2_8175080", 9) 
squash <- movemarker(squash, "2_8224962", 9) 
plotRF(squash, chr = c(2,9))

#########################          Estimate final map order and distances        ################################

#Now plot recombination fractions again
pdf("plots/RF_correct_LGs.pdf")
plotRF(squash)
dev.off()

###### Use ASMap to reorder markers within each linkage group
squash.mst <- convert2bcsft(squash, BC.gen = 0, F.gen = 0, error.prob = 0.003)
squash.mst <- mstmap.cross(squash.mst, id='Sample', bychr=T, p.value=2, dist.fun="kosambi")

#Invert all marker orders so genetic and physical maps are in same orientation
squash.mst <- flip.order(squash.mst, 1:20)

#Jitter marker positions
squash.mst <- jittermap(squash.mst)

pdf("plots/genetic_map.pdf")
plot.map(squash.mst)
dev.off()

## Plot physical distance vs genetic distance
pdf("plots/gen_vs_phys.pdf", width=6, height=4)
old.par <- par(no.readonly = T)
par(mfrow=c(4,5), mar=c(2,2,1,1), oma=c(4,4,1,1), xpd=NA)
for(i in 1:20){
  xlab <- ''
  ylab <- ''
  if(i == 18){
    xlab <- "    Physical position (Mb)"
  }
  if(i == 11){
    ylab <- "                   Genetic position (cM)"
  }
  cm <- pull.map(squash.mst, i)
  cm <- unlist(cm)
  bp <- as.integer(gsub(paste(i, '.' , i, '_',sep=''), '', names(cm)))
  plot(bp, cm, xaxt='n', yaxt='n',
       xlim=c(0,max(bp, na.rm=T)),
       xlab=xlab, 
       ylab=ylab,
       main = paste("LG",i),
       cex.main=0.9, cex.lab=1.1)
  axis(1, 
       at=seq(0,max(bp, na.rm=T),by=1000000),
       labels=seq(0,max(bp, na.rm=T),by=1000000)/1000000,
       line=0, cex.axis=0.75, padj=0, tcl=-0.3, mgp=c(3,0.25,0))
  axis(2, 
       at=seq(0,max(cm),by=15),
       labels=seq(0,max(cm),by=15), las=2,
       line=0, cex.axis=0.75, padj=0, tcl=-0.3, mgp=c(3,0.50,0))
}
par(old.par)
dev.off()

squash.mst <- calc.genoprob(squash, step=1, error.prob = .003, map.function = "kosambi", stepwidth = "fixed")

saveRDS(squash.mst, "squash_map.rds")










