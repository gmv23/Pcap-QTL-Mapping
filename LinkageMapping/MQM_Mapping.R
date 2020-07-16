setwd("~/Documents/work/Smart_lab/P_capsici/xqtl/F23_validation/analysis/LinkageMapping/")
library(qtl)

squash <- readRDS("squash_map.rds")
squash <- calc.genoprob(squash, error.prob = .003)

perms.out <- readRDS(paste("HK_permutations/perms1.rds"))
for(i in c(2:10)){ 
  perms.i <- readRDS(paste("HK_permutations/perms", i, ".rds", sep=""))
  print(i)
  print(summary(perms.i))
  perms.out <- c(perms.out, perms.i)
}
summary(perms.out)

pen <- calc.penalties(perms.out)
stepwise.out <- stepwiseqtl(squash, max.qtl=8, penalties=pen, verbose=T, 
                               pheno.col = 5, method = "hk", 
                               keeplodprofile = T, keeptrace = T)

qtl.fit <- fitqtl(squash, qtl=stepwise.out, get.ests = T, pheno.col=5, method='hk')
summary(qtl.fit)

#Pull out imputed genotypes
geno <- pull.argmaxgeno(argmax.geno(squash))
geno <- geno - 1
rownames(geno) <- squash$pheno$Sample
pheno <- squash$pheno
write.csv(pheno, "pheno.csv", quote=F, row.names = F)


#Prune markers
markers_to_remove <- c()
gen_map <- pull.map(squash)
for(lg in 1:20){
  map.lg <- gen_map[[lg]]
  markers <- names(map.lg)
  last_distance <- map.lg[1]
  for(i in 2:length(map.lg)){
    if(map.lg[i] - last_distance < 1){
      markers_to_remove <- c(markers_to_remove, markers[i])
    }
      last_distance <- map.lg[i]
  }
}

geno_prune <- geno[,!colnames(geno) %in% markers_to_remove]
write.csv(geno, "geno_imp.csv", row.names = T)
write.csv(geno_prune, "geno_prune.csv", row.names = T)


snps <- data.frame("CHR" = sapply(colnames(geno), function(x) unlist(strsplit(x,"_"))[1]),
                   "BP" = as.integer(sapply(colnames(geno), function(x) unlist(strsplit(x,"_"))[2])))

#Bring in multipool results and get position of peak on Chroms 4,5,8,16
lods <- read.table("../../../analysis/round2/multipool_results/joint_multi_lods.txt")
colnames(lods) <- c("CHROM", "POS", "Rep1", "Rep2", "Pool")
lods$CHROM <- as.integer(gsub("Cp4.1LG", "",lods$CHROM))
peaks <- aggregate(cbind(Rep1, Rep2) ~ CHROM, data=lods, FUN=max)
#Now get closest GBS marker to peak
marker_pos <- rep(NA, 4)
lod_peaks <- rep(NA, 4)
chroms <- c(4,5,8,12,16)
for(i in 1:length(chroms)){
  chrom <- chroms[i]
  peak_lod1 <- mean(lods$POS[which(lods$Rep1==peaks$Rep1[peaks$CHROM==chrom] & lods$CHROM==chrom)])
  peak_lod2 <- mean(lods$POS[which(lods$Rep2==peaks$Rep2[peaks$CHROM==chrom] & lods$CHROM==chrom)])
  peak_lod <- mean(c(peak_lod1, peak_lod2)) #Mean of both LODs
  lod_peaks[i] <- peak_lod
  chrom_meta <- snps[snps$CHR==chrom,]
  chrom_diff <- abs(chrom_meta$BP-peak_lod)
  peak_gbs <- chrom_meta[which(chrom_diff==min(chrom_diff)),]
  marker_pos[i] <- which(snps$CHR==chrom & snps$BP==peak_gbs$BP)
}
marker_distances <- abs(snps$BP[marker_pos] - lod_peaks); marker_distances #Distance from GBS marker to LOD peak

qtl_geno <- geno[,marker_pos]
write.csv(qtl_geno, 'qtl_geno.csv', row.names = T)


disease <- data.frame("Sample" = squash$pheno$Sample, 
                 "audpc" = squash$pheno$audpc, 
                 geno[,marker_pos])
rownames(disease) <- NULL

disease_orient <- disease
disease_orient$X16_7011538 <- abs(disease_orient$X16_7011538 - 2)
boxplot(disease_orient$audpc ~ disease_orient$X4_9371608)
boxplot(disease_orient$audpc ~ disease_orient$X5_682457)
boxplot(disease_orient$audpc ~ disease_orient$X8_2334412)

number_markers <- apply(disease_orient[,3:6], 1, function(x) sum(x==2, na.rm=T))
pdf("plots/homoQTL.pdf")
boxplot(disease$audpc ~ number_markers,
        xlab = "Number of markers homozygous for beneficial allele",
        ylab = "AUDPC BLUP")
dev.off()

#Allele frequencies in random vs selected
random <- 1:161
selected <- 162:180
af <- function(x){sum(x)/(2*length(x))}
random_afs <- apply(disease[random,3:6],2,af)
selected_afs <- apply(disease[selected,3:6],2,af)

pdf("plots/qtl_on_map.pdf", height = 6, width = 4)
plot(stepwise.out, chr=c(4,5,8,19))
dev.off()

bsa_marker_pos <- find.markerpos(squash, c('4_9371608', '5_682457', '8_2334412', '12_2429191', '16_7011538'))
bsa_qtl <- makeqtl(squash, chr=c(4,5,8,12,16), pos=bsa_marker_pos$pos, what = "prob")
lm_qtl <- stepwise.out
bsa_qtl.fit <- fitqtl(squash, pheno.col = 5, qtl=bsa_qtl, 
                      formula = Y~Q1+Q2+Q3+Q4+Q5, method = 'hk')
lm_qtl.fit <- fitqtl(squash, pheno.col = 5, qtl= stepwise.out, 
                     formula = Y~Q1+Q2+Q3+Q4, method = 'hk')
summary(bsa_qtl.fit)
summary(lm_qtl.fit)

#Look at GP vs QTL

#Get rid of samples with no phenotype
no_pheno <- disease$Sample[is.na(disease$audpc)]
disease <- disease[!disease$Sample %in% no_pheno,]
geno <- geno[!rownames(geno) %in% no_pheno, !colnames(geno) %in% no_pheno]

#Fit genomewide prediction model
library(rrBLUP)
K <- A.mat(geno-1)
pred <- kin.blup(data = disease, geno = "Sample", pheno = "audpc", K=K)
pred$Vg/(pred$Vg+pred$Ve)
n <- nrow(disease)

#Compare GP to OLS
accuracies_ols <- rep(NA, 100)
accuracies_gblup <- rep(NA, 100)
set.seed(801)
for(i in 1:100){
  sets <- rep(1:5, n)[sample(1:n, size=n, replace=F)]
  fivefold_ols <- rep(NA,5)
  fivefold_gblup <- rep(NA,5)
  for(j in 1:5){
    #GBLUP
    masked <- which(sets==j)
    pheno <- disease
    pheno[masked,'audpc'] <- NA
    rrblup.out <- kin.blup(data = pheno, geno = "Sample", pheno = "audpc", K=K)
    gebvs <- rrblup.out$g
    accuracy_gblup <- cor(disease$audpc[masked], gebvs[masked])
    fivefold_gblup[j] <- accuracy_gblup
    
    #OLS
    pheno.train <- disease[sets!=j,]
    pheno.test <- disease[sets==j,]
    train.lm <- lm(audpc ~ X4_9371608 + X5_682457 + X8_2334412 + X16_7011538, data=pheno.train)
    blups.pred <- predict(train.lm, pheno.test)
    blups.test <- pheno.test$audpc
    accuracy_ols <- round(cor(blups.pred, blups.test),2)
    fivefold_ols[j] <- accuracy_ols
  }
  accuracies_gblup[i] <- mean(fivefold_gblup)
  accuracies_ols[i] <- mean(fivefold_ols)
  
}
pdf("plots/pred_accuracies.pdf")
boxplot(accuracies_ols, accuracies_gblup, 
        ylab = 'cor(predicted, observed)',
        names = c('OLS', 'GBLUP'))
dev.off()

#### Rate families based on markers, phenos, or gebvs
rrblup.out <- kin.blup(data = pheno, geno = "Sample", pheno = "audpc", K=K)
gebvs <- sort(rrblup.out$g)

phenos <- disease$audpc
names(phenos) <- disease$Sample
phenos <- sort(phenos)

marker_scores <- apply(disease_orient[,3:6], 1, sum)
names(marker_scores) <- disease_orient$Sample
marker_scores <- sort(marker_scores, decreasing =T)

markers.lm <- lm(audpc ~ X4_9371608 + X5_682457 + X8_2334412 + X16_7011538, data=disease)
marker_predictions <- predict(markers.lm)
names(marker_predictions) <- disease$Sample
marker_predictions <- sort(marker_predictions)

best_families <- cbind(names(head(phenos,10)), 
      names(head(gebvs,10)), 
      names(head(marker_scores,10)), 
      names(head(marker_predictions,10)))

write.csv(best_families, "best_families.csv")
