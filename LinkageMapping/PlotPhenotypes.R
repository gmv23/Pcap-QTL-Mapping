setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/LinkageMapping/")

####################################     Import and clean data    #################################

blues <- read.csv("../F23/tables/f23_blues.csv", stringsAsFactors = F)
colnames(blues)[1] <- "Sample"

geno <- read.csv("tables/qtl_marker_genotypes.csv", stringsAsFactors = F)
colnames(geno)[1] <- "Sample"

####################################     Plot phenotypic distributions    #################################

pdf("plots/phenotypic_distributions.pdf", width=7, height = 6)

#Get plot layout
m <- rbind(c(1,1,2,3,3,4),
           c(5))
layout(m)

#Set par settings
old.par <- par(no.readonly = T)
par(oma=c(1,1,1,1), 
    mar=c(5,4,2,1)
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
text(x=-35,y=-20,labels="C", font=1,cex=2)
text(x=-35,y=45,labels="A", font=1,cex=2)
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
text(x=-35,y=45,labels="B", font=1,cex=2)

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
#Plot 3: phenotypes depending on QTL genotype

#Get rid of checks
blues.f2 <- blues[blues$Sample %in% geno$Sample,]
geno <- geno[match(blues.f2$Sample, geno$Sample),]

#Get list with phenotypes for each QTL allele
blues.cond <- list()
for(i in 1:5){
  for(j in 0:2){
    blues.ij <- blues.f2$raudpc[which(geno[,(i+1)] == j)]
    blues.name <- paste(i,j,sep="_")
    blues.cond[[blues.name]] <- blues.ij
  }
}

#Conduct Tukeys HSD for each QTL comparing alleles
letter_assignments <- c()
for(i in 1:5){
  data = data.frame("raudpc" = blues.f2$raudpc,
                    "geno" = as.factor(geno[,(i+1)]))
  data.lm <- lm(raudpc~geno,data=data)
  data.aov <- aov(data.lm)
  data.hsd <- HSD.test(data.aov, "geno")
  data.groups <- data.hsd$groups
  data.groups <- data.groups[order(row.names(data.groups)),]
  letter_assignments <- c(letter_assignments, as.character(data.groups$groups))
}

#Where will the boxplots be
at <- c(1,2,3,6,7,8,11,12,13,16,17,18,21,22,23)

#Get boxplot coordinates
boxplot.data <- boxplot(blues.cond, plot=F)

#Draw boxplot
boxplot(blues.cond,at=at,xaxt='n', 
        col=rep(brewer.pal(3,"Dark2"),5),
        ylab = "rAUDPC",
        cex.lab=1.1,
        ylim=c(min(blues.f2$raudpc),100))

#Add Tukey letters
text(x=at,y=boxplot.data$stats[5,]+5,labels=letter_assignments)

#Add QTL names to x axis
axis(side=1,at=c(2,7,12,17,22),labels=
       gsub("X","",colnames(geno)[-1]),
     cex.axis=1.2)

#Plot 3 legend
par(xpd=NA)
legend(1,-1, bty='n', ncol=3,
       fill = brewer.pal(3,"Dark2"),
       legend = c("Homozygous Dunja",
                  "Heterozygous",
                  "Homozygous Pc-NY21"),
       cex=1.1,
       x.intersp=0.6,text.width=6.5)

par(old.par)
dev.off()
