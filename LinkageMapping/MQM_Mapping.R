setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/LinkageMapping/")
library(qtl)

#Import squash genetic map 
squash <- readRDS("../genetic_map/squash_map.rds")

#Import snps
snps <- read.csv("../genetic_map/tables/snps_filter_prune.csv")

########################################      Perform MQM using F23 data    ####################################

#Import LOD scores from 1000 
#permutations of 2-d genomewide scan to get model penalties
perms.out <- readRDS(paste("data/out/perms1.rds"))
for(i in c(2:50)){ 
  perms.i <- readRDS(paste("data/out/perms", i, ".rds", sep=""))
  print(i)
  print(summary(perms.i))
  perms.out <- c(perms.out, perms.i)
}
summary(perms.out)

#Calculate penalties for model selection
pen <- calc.penalties(perms.out)

#Perform automatic model selection
stepwise.out <- stepwiseqtl(squash, penalties=pen, verbose=T, 
                               pheno.col = 5, method = "hk", 
                               keeplodprofile = T, keeptrace = T)

#Plot LOD profile of QTL
pdf("plots/MQM_LOD_profile.pdf", width=7, height=5)
plotLodProfile(stepwise.out, ylab = "LOD")
dev.off()

#Get 90% confidence intervals for QTL
CI.lm <- data.frame("Chrom" = stepwise.out$chr,
                      "Start" = NA,
                      "Peak" = NA,
                      "End" = NA)
for(i in 1:nrow(CI.lm)){
   ci <- bayesint(stepwise.out, qtl.index=i, expandtomarkers = T, prob=0.90)
   marker.pos <- find.marker(squash, chr=ci$chr, pos=ci$pos)
   bp.pos <- as.integer(gsub("[0-9]?_","", marker.pos))
   CI.lm[i,2:ncol(CI.lm)] <- bp.pos  
}
  

###################################      Find GBS markers closest to BSA peaks    #################################

#Pull out imputed genotypes
geno <- pull.argmaxgeno(argmax.geno(squash))
geno <- geno - 1
rownames(geno) <- squash$pheno$Sample

#Bring in multipool results to find GBS marker position closest to BSA Lod peaks on Chroms 4,5,8,12,16
peaks <- read.csv("../BSA/data/CIs_and_peaks.txt", header=T)

#Make new data frame to put results with GBS markers
BSA_markers <- peaks[c(4,5,8,12,16),c("CHROM", "Rep1_peak", "Rep2_peak")]
BSA_markers$meanpeak <- apply(cbind(BSA_markers$Rep1_peak, BSA_markers$Rep2_peak),1,mean)
BSA_markers$GBS_BP <- NA

#Now get closest GBS marker to peak
for(i in 1:nrow(BSA_markers)){
  chrom <- BSA_markers$CHROM[i]
  chrom_meta <- snps[snps$CHROM==chrom,]
  chrom_diff <- abs(chrom_meta$BP-BSA_markers$meanpeak[i])
  peak_gbs <- chrom_meta[which(chrom_diff==min(chrom_diff)),]
  BSA_markers$GBS_BP[i] <- peak_gbs$BP
}

BSA_markers$distance <- abs(BSA_markers$meanpeak - BSA_markers$GBS_BP)
BSA_markers$marker_name <- paste(as.integer(gsub("Cp4.1LG", "", BSA_markers$CHROM)),
                                 BSA_markers$GBS_BP, sep = "_")

#Pull out genotypes at the marker loci
qtl_geno <- geno[,BSA_markers$marker_name]

#Look at Pc-NY21 allele frequency in 'random' and 'selected' F2 individuals
which.random <- grep("Vog",rownames(geno))
which.select <- grep("18_G", rownames(geno))
BSA_markers$AF <- apply(qtl_geno, 2, function(x) sum(x/(length(x)*2)))
BSA_markers$AF_random <- apply(qtl_geno[which.random,], 2, function(x) sum(x/(length(x)*2)))
BSA_markers$AF_select <- apply(qtl_geno[which.select,], 2, function(x) sum(x/(length(x)*2)))

#Write table
GBS_markers_print <- BSA_markers[,c("CHROM", "GBS_BP", "distance", "AF_random", "AF_select")]
write.csv(GBS_markers_print, "tables/GBS_markers_print", row.names = F, quote = F)

####################################      Look at QTL effects for two models    #################################

#Linkage mapping QTL model
qtl.lm <- stepwise.out

#BSA QTL model
bsa_marker_pos <- find.markerpos(squash, BSA_markers$marker_name)
qtl.bsa <- makeqtl(squash, chr = bsa_marker_pos$chr, pos = bsa_marker_pos$pos, what="prob")

#Fit QTL models and look at % variation and effect sizes
fit.lm <- fitqtl(squash, pheno.col=5, qtl=qtl.lm, method="hk")
fit.bsa <- fitqtl(squash, pheno.col=5, qtl=qtl.bsa, method="hk")











####################################      Do BSA in individually genotyped F2s    #################################

#chisquare test for allele frequency differences between random and select individuals

count_alleles <- function(x){
  n1 <- sum(x==1)+2*sum(x==0)
  n2 <- sum(x==1)+2*sum(x==2)
  return(c(n1,n2))
}
pvals <- rep(NA, ncol(geno))
for(i in 1:ncol(geno)){
  alleles.ran <- count_alleles(geno[which.random,i])
  alleles.sel <- count_alleles(geno[which.select,i])
  con.table <- cbind(alleles.ran, alleles.sel)
  chisq <- chisq.test(con.table) 
  pvals[i] <- chisq$p.value
}

chisq.results <- data.frame("CHR"=as.integer(gsub("Cp4.1LG", "", snps$CHROM)),
                            "BP"=snps$BP,
                            "P" = pvals)














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