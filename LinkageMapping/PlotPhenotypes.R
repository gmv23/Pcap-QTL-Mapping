setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/LinkageMapping/")

####################################     Import and clean data    #################################

blues <- read.csv("../F23/tables/f23_blues.csv", stringsAsFactors = F)
colnames(blues)[1] <- "Sample"

geno <- read.csv("tables/qtl_marker_genotypes.csv", stringsAsFactors = F)
colnames(geno)[1] <- "Sample"

####################################     Plot phenotypic distributions    #################################

pdf("plots/phenotypic_distributions.pdf", width=7, height = 6)

#Get plot layout
m <- rbind(rep(c(1,1,2,3,3,4),each=5),
           rep(c(5,6,7,8,9),each=6))
layout(m)

#Set par settings
old.par <- par(no.readonly = T)
par(oma=c(3,3,1,1), 
    mar=c(6,2,3,3)
)

#Get data subsets for histogram
which.ran <- grep("Vog", blues$Sample)
which.sel <- grep("18_G", blues$Sample)
which.dunja <- grep("Dunja", blues$Sample)
which.pcny21 <- grep("17_G_14", blues$Sample)
which.f1 <- grep("17_G_020", blues$Sample)

#Plot 1: phenotypic distribution and checks

hist(blues$raudpc[c(which.ran, which.sel)],
     xlim = c(0,100),
     ylim=c(0,40),
     breaks=seq(0,100,by=5),
     col=adjustcolor(colors()[543],alpha.f=0.75),
     xlab = "rAUDPC",
     ylab = "Frequency",
     main = "",
     border=colors()[185],
     cex.lab=1.1)
abline(v=blues$raudpc[which.dunja], lwd=2, lty = 1)
abline(v=blues$raudpc[which.f1], lwd=2, lty=2)
abline(v=blues$raudpc[which.pcny21], lwd=2, lty=3)

#Plot 1 legend
par(xpd=NA)
text(x=-35,y=45,labels="A", font=1,cex=1.5)
plot(0, type='n', xaxt='n', yaxt='n',
     xlab='', ylab='', bty='n',
     xlim=c(0,10),ylim=c(0,10))
legend(-17,10, bty='n',
       lty=c(1,2,3), lwd=2,
       legend = c("Dunja", 
                  expression(atop("Pc-NY21",
                                  paste("x Dunja ", "F"["1"]))), 
                  "Pc-NY21"),
       y.intersp=2)

par(xpd=F)
#Plot 2: phenotypic distribution of random vs selected F23s

hist(blues$raudpc[c(which.ran)],
     xlim = c(0,100),
     ylim=c(0,40),
     breaks=seq(0,100,by=5),
     col=adjustcolor(colors()[430],alpha.f=0.75),
     xlab = "rAUDPC",
     ylab = "Frequency",
     main = "",
     border=colors()[185],
     cex.lab=1.1)
hist(blues$raudpc[c(which.sel)],
     breaks=seq(0,100,by=5),
     col=adjustcolor(colors()[630],alpha.f=0.75),
     main = "", 
     border=colors()[185],
     add=T)

points(x=rep(median(blues$raudpc[which.ran]),25), y = seq(0,40,length.out = 25),
       cex=0.7,bg=colors()[430], col=colors()[190], pch=22)
points(x=rep(median(blues$raudpc[which.sel]),25), y = seq(0,40,length.out = 25),
       cex=0.7,bg=colors()[630], col=colors()[190], pch=22)
#Plot 2 legend
par(xpd=NA)
text(x=-35,y=45,labels="B", font=1,cex=1.5)

plot(0, type='n', xaxt='n', yaxt='n',
     xlab='', ylab='', bty='n',
     xlim=c(0,10),ylim=c(0,10))
legend(-17,10, bty='n',
       pch=22, pt.bg=c(colors()[430],colors()[630]), col=colors()[190],
       pt.cex=2,
       legend = c(expression(paste("Random ", "F"["2:3"])),
                  expression(paste("Selected ", "F"["2:3"]))),
       y.intersp=2)
legend(-20,7, bty='n',
       pch=22, pt.bg=c(colors()[430],colors()[630]), col=colors()[190],
       legend = c(expression(atop(" Median",
                                  paste(" random ", "F"["2:3"]))),
                  expression(atop(" Median",
                                  paste(" selected ", "F"["2:3"])))),
       y.intersp=2.5, text.col='white')
legend(-18,7, bty='n',
       pch=22, pt.bg=c(colors()[430],colors()[630]), col=colors()[190],
       legend = c(expression(atop(" Median",
                                  paste(" random ", "F"["2:3"]))),
                  expression(atop(" Median",
                                  paste(" selected ", "F"["2:3"])))),
       y.intersp=2.5, text.col='white')
legend(-16,7, bty='n',
       pch=22, pt.bg=c(colors()[430],colors()[630]), col=colors()[190],
       legend = c(expression(atop(" Median",
                                  paste(" random ", "F"["2:3"]))),
                  expression(atop("  Median",
                                  paste(" selected ", "F"["2:3"])))),
       y.intersp=2.5)
par(xpd=F)

#Plot 3-7: phenotypes depending on QTL genotype

#New par
par(mar=c(5,1.75,4,1.75), bty='n', xpd=NA)

#Get rid of checks
blues.f2 <- blues[blues$Sample %in% geno$Sample,]
geno <- geno[match(blues.f2$Sample, geno$Sample),]

for(i in 1:5){
  qtl.geno <- geno[,(i+1)]
  #Get boxplot coordinates
  boxplot.data <- boxplot(blues.f2$raudpc~qtl.geno, plot=F)
  
  if(i == 1){
    ylab <- "rAUDPC"
    xlab <- ""
  }else if(i == 3){
    xlab <- "Marker genotype"
    ylab <- ""
  }else{
    ylab <- ""
    xlab <- ""
  }
  #Draw boxplot
  boxplot(blues.f2$raudpc~qtl.geno,
          names=c("AA","AB","BB"),
          col=brewer.pal(3,"Dark2"),
          ylab = ylab, xlab=xlab,
          cex.lab=1.1,
          ylim=c(min(blues.f2$raudpc),100),
          yaxt = 's')
  
  #Add marker name
  text(2,110,
       label=gsub("X","",colnames(geno)[i+1]),
       font=4)
  
  #Add letter
  text(-0.6,120,LETTERS[(i+2)], cex=1.5)
  
  #Tukeys
  data.lm <- lm(blues.f2$raudpc~qtl.geno)
  data.aov <- aov(data.lm)
  data.hsd <- HSD.test(data.aov, "qtl.geno")
  data.groups <- data.hsd$groups
  data.groups <- data.groups[order(row.names(data.groups)),]
  
  #Add Tukey letters
  text(x=1:3,
       y=boxplot.data$stats[5,]+5,
       labels=as.character(data.groups$groups))
  
  if(i==3){
    legend(-9,-18, bty='n', ncol=3,
           fill = brewer.pal(3,"Dark2"),
           legend = c("AA = Homozygous Dunja",
                      "AB = Heterozygous",
                      "BB = Homozygous Pc-NY21"),
           cex=1.1,
           x.intersp=0.6,text.width=6.5)
  }
  
}

par(old.par)
dev.off()
