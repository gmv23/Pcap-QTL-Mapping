setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/PredictionModels/")
library(rrBLUP)
library(BGLR)
library(RColorBrewer)

##################################   Clean and load data   #############################

geno <- read.csv("../LinkageMapping/tables/imputed_genos.csv", row.names = 1)
pheno <- read.csv("../genetic_map/tables/blues_final.csv")
qtl_geno <- read.csv("../LinkageMapping/tables/qtl_marker_genotypes.csv", row.names=1)


#Get rid of NAs and get in same order
pheno <- pheno[-which(is.na(pheno$raudpc)),]
sample_intersection <- intersect(pheno$Sample, rownames(geno))

geno <- geno[match(sample_intersection, rownames(geno)),]
pheno <- pheno[match(sample_intersection, pheno$Sample),]
qtl_geno <- qtl_geno[match(sample_intersection, rownames(qtl_geno)),]

#Add AUDPC to qtl_geno object
qtl_geno$raudpc <- pheno$raudpc

#Get genomic relationship matrix
K <- A.mat(geno-1)

#Get narrow sense heritability
pred <- kin.blup(data = pheno, geno = "Sample", pheno = "raudpc", K=K)
pred$Vg/(pred$Vg+pred$Ve)

##################################   Run cross-validation   #############################

#Compare 4 models:
# 1) LS using 4 markers in each of 4 BSA QTL
# 2) GBLUP
# 3) GBLUP with marker fixed effects
# 4) Bayes-B

accuracies <- matrix(NA, nrow=50, ncol=4)
n <- nrow(geno)
set.seed(801)
counter <- 0
for(i in 1:50){ #20 replications
  print(paste('Replication', i))
  #Get folds (fivefold)
  which.train <- sample(n, round(0.8*n), replace=F)
  which.test <- (1:n)[!(1:n) %in% which.train]

  #Least squares
  pheno.train <- qtl_geno[which.train,]
  pheno.test <- qtl_geno[which.test,]
  train.lm <- lm(raudpc ~ 1 + X4_8908173 + X5_385710 + X8_3462310 + X16_7011538 + X12_2386018, data=pheno.train)
  blues.pred <- predict(train.lm, pheno.test)
  blues.test <- pheno.test$raudpc
  accuracy <- cor(blues.pred, blues.test)
  accuracies[i,1] <- accuracy
    
  #GBLUP
  pheno.fold <- pheno[,c("Sample", "raudpc")]
  pheno.fold[which.test,'raudpc'] <- NA
  pheno.fold$intercept <- 1
  rrblup.out <- kin.blup(data = pheno.fold, geno = "Sample", pheno = "raudpc", K=K, fixed="intercept")
  gebvs <- rrblup.out$g
  accuracy <- cor(pheno$raudpc[which.test], gebvs[which.test])
  accuracies[i,2] <- accuracy
    
  #GBLUP with fixed marker effects
  Z <- diag(nrow=nrow(pheno))
  X <- as.matrix(cbind(1,qtl_geno[,1:5]))
  rrblup.out <- mixed.solve(y = pheno.fold$raudpc, Z = Z,
                            K=K, X = X, method="REML")
  gebvs <- X %*% rrblup.out$beta + Z %*% rrblup.out$u
  accuracy <- cor(pheno$raudpc[which.test], gebvs[which.test])
  accuracies[i,3] <- accuracy
    
  #Bayes-B
  bglr.out <- BGLR(y=pheno.fold$raudpc, 
                   ETA=list(list(X=geno,model="BayesB")), verbose=F)
  gebvs <- bglr.out$yHat
  accuracy <- cor(pheno$raudpc[which.test], gebvs[which.test])
  accuracies[i,4] <- accuracy
}

##################################  Make plot   #############################

pdf("plots/accuracies.pdf", width=3.30709, height=4)
old.par <- par(no.readonly = T)
par(mgp=c(3,1,0), oma=c(7,4,2,0.5), mar=c(0,0,0,0), xpd=NA)
plot(0, type='n',xaxt='n',yaxt='n',ylim=c(0,max(accuracies)), xlim=c(0.6,4.4),
     xlab="", ylab="Pearson's cor(predicted,observed)")
par(xpd=T)
for(i in seq(0,1,by=0.1)){
  abline(h=i, col='gray', lty=2)
}
par(xpd=NA)
boxplot(accuracies, las=2,
        names = c("QTL MLR", "GBLUP", "GBLUP + QTL", "Bayes-B"),
        col = brewer.pal(4, "Set2"),
        cex.lab=1, add=T)
#        names=rep("",4))

#axis(1, at=1:4, 
#     labels = c("QTL\nMLR", "GBLUP", "GBLUP\n + QTL", "Bayes-B"), las=2,
#     cex.axis=0.6, padj=0)
par(old.par)
dev.off()



apply(accuracies,2,median)

