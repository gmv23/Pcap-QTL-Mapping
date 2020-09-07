setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/LinkageMapping/")
library(qtl)

#Import squash genetic map 
squash <- readRDS("../genetic_map/squash_map.rds")
#squash$pheno$raudpc[163:181] <- NA #to look at results only using 'random' f23s

#Import snps
snps <- read.csv("../genetic_map/tables/snps_filter_prune.csv")

#Import multipool results to find GBS markers closest to BSA peaks
peaks <- read.csv("../BSA/data/CIs_and_peaks.txt", header=T)

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
pdf("plots/MQM_LOD_profile.pdf", width=7, height=7)
jpeg("plots/MQM_LOD_profile.jpeg", width=7, height=7, units="in", res=600)
#Rename QTL for purpose of plot
stepwise.out.rename <- stepwise.out
stepwise.out.rename$name <- paste(stepwise.out.rename$name, "cM", sep=" ")
names(attr(stepwise.out.rename, "lodprofile")) <- paste(names(attr(stepwise.out.rename, "lodprofile")), "cM", sep=" ")
plotLodProfile(stepwise.out.rename, ylab = "LOD")
dev.off()

#Get 90% confidence intervals for QTL
CI.lm <- data.frame("Chrom" = stepwise.out$chr,
                      "Peak" = NA,
                      "CI" = NA)
for(i in 1:nrow(CI.lm)){
   ci <- bayesint(stepwise.out, qtl.index=i, expandtomarkers = T, prob=0.90)
   marker.pos <- find.marker(squash, chr=ci$chr, pos=ci$pos)
   bp.pos <- as.integer(gsub("[0-9]+_","", marker.pos))
   ci.string <- paste(bp.pos[1],bp.pos[3],sep="-")
   CI.lm$Peak[i] <- bp.pos[2]
   CI.lm$CI[i] <- ci.string
}
#Hardcode in Chrom 4 CI; disjoint interval because of inversion in genetic map
CI.lm$CI[1] <- "0-20944;8074112-8504971"

#For comparison results --- look at LOD scores and CIs from 1D QTL scan
#out.perm <- scanone(squash, n.perm=1000, pheno.col=5, verbose=FALSE)
summary(out.perm)
out.1d <- scanone(squash, pheno.col=5, method = "hk")
plot(out.1d)

CI.1d <- data.frame("Chrom" = stepwise.out$chr,
                    "Peak" = NA,
                    "CI" = NA)
for(i in 1:nrow(CI.1d)){
  ci <- bayesint(out.1d, chr=as.integer(as.character(CI.1d$Chrom))[i], expandtomarkers = T, prob=0.90)
  marker.pos <- find.marker(squash, chr=ci$chr, pos=ci$pos)
  bp.pos <- as.integer(gsub("[0-9]+_","", marker.pos))
  ci.string <- paste(bp.pos[1],bp.pos[3],sep="-")
  CI.1d$Peak[i] <- bp.pos[2]
  CI.1d$CI[i] <- ci.string
}
###################################      Find GBS markers closest to BSA peaks    #################################

#Pull out imputed genotypes
geno <- pull.argmaxgeno(argmax.geno(squash, error.prob=0.001,
                        map.function="kosambi", stepwidth="fixed"))
geno <- geno - 1
rownames(geno) <- squash$pheno$Sample

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

#Write tables
GBS_markers_print <- BSA_markers[,c("CHROM", "GBS_BP", "distance", "AF_random", "AF_select")]
GBS_markers_print$CHROM <- as.integer(gsub("Cp4.1LG","",as.character(GBS_markers_print$CHROM)))
write.csv(GBS_markers_print, "tables/GBS_markers_print.csv", row.names = F, quote = F)
write.csv(qtl_geno, "tables/qtl_marker_genotypes.csv", row.names = T, quote=F)
write.csv(geno, "tables/imputed_genos.csv", row.names=T, quote=F)

####################################      Look at QTL effects for two models    #################################

############## Linkage mapping QTL model ###########
qtl.lm <- stepwise.out

#Fit QTL models and look at % variation and effect sizes
fit.lm <- fitqtl(squash, pheno.col=5, qtl=qtl.lm, method="hk", get.ests=T)

#Total % variation explained by QTL model
fit.lm$result.full

#QTl effect sizes
effects.lm <- fit.lm$ests$ests

#Put together QTL location, CI, % variation explained, p-value, and effect sizes
qtl_summary.lm <- cbind(CI.lm,
                        "Genetic_pos" = qtl.lm$pos,
                        fit.lm$result.drop[,c("LOD", "%var", "Pvalue(F)")])
qtl_summary.lm <- data.frame(qtl_summary.lm)
#Add effect sizes
qtl_summary.lm$effect_a <- NA
qtl_summary.lm$effect_d <- NA
for(i in 1:nrow(qtl_summary.lm)){
  qtl_summary.lm[i,c("effect_a", "effect_d")] <- 
    effects.lm[grep(paste(qtl_summary.lm$Chrom[i],"@",sep=""), names(effects.lm))]
}

############## BSA QTL model ###########

bsa_marker_pos <- find.markerpos(squash, BSA_markers$marker_name)
qtl.bsa <- makeqtl(squash, chr = bsa_marker_pos$chr, pos = bsa_marker_pos$pos, what="prob")

#Fit QTL models and look at % variation and effect sizes
fit.bsa <- fitqtl(squash, pheno.col=5, qtl=qtl.bsa, method="hk", get.ests = T)

#Total % variation explained by QTL model
fit.bsa$result.full

#QTl effect sizes
effects.bsa <- fit.bsa$ests$ests

#Put together QTL location, CI, % variation explained, p-value, and effect sizes
qtl_summary.bsa <- cbind("Chrom" = GBS_markers_print$CHROM,
                        "Position" = GBS_markers_print$GBS_BP,
                        "Genetic_pos" = qtl.bsa$pos,
                        fit.bsa$result.drop[,c("LOD", "%var", "Pvalue(F)")])
qtl_summary.bsa <- data.frame(qtl_summary.bsa)
#Add effect sizes
qtl_summary.bsa$effect_a <- NA
qtl_summary.bsa$effect_d <- NA
for(i in 1:nrow(qtl_summary.bsa)){
  qtl_summary.bsa[i,c("effect_a", "effect_d")] <- 
    effects.bsa[grep(paste(qtl_summary.bsa$Chrom[i],"@",sep=""), names(effects.bsa))]
}

####################################      Write QTL tables    #################################

qtl_summary.bsa.round <- qtl_summary.bsa
qtl_summary.bsa.round[,c(3,4,5,7,8)] <- apply(qtl_summary.bsa[,c(3,4,5,7,8)],2,round,digits=2)
qtl_summary.bsa.round$Pvalue.F. <- signif(qtl_summary.bsa.round$Pvalue.F.,3)
qtl_summary.bsa.round$Position <- round(qtl_summary.bsa.round$Position/1000000,2)

qtl_summary.lm.round <- qtl_summary.lm
qtl_summary.lm.round[,c(4,5,6,8,9)] <- apply(qtl_summary.lm[,c(4,5,6,8,9)],2,round,digits=2)
qtl_summary.lm.round$Pvalue.F. <- signif(qtl_summary.lm.round$Pvalue.F.,3)

#Turn numbers in CI into Mb
round_ci <- function(x){
  
  m.numbers <- gregexpr('[0-9]+',x)
  numbers <- as.integer(unlist(regmatches(x,m.numbers)))
  numbers.round <- round(numbers/1000000,2)
  
  m.connectors <- gregexpr('[^0-9]+',x)
  connectors <- unlist(regmatches(x,m.connectors))
  
  newstring <- ""
  for(i in 1:length(connectors)){
    newstring <- paste(newstring,numbers.round[i], sep="")
    newstring <- paste(newstring,connectors[i], sep="")
    if(i == length(connectors)){
      newstring <- paste(newstring,numbers.round[(i+1)], sep="")
    }
  }
  return(newstring)
}

qtl_summary.lm.round$CI <- sapply(qtl_summary.lm.round$CI, round_ci)
qtl_summary.lm.round$Peak <- round(qtl_summary.lm.round$Peak/1000000,2)

write.csv(qtl_summary.bsa, "tables/qtl_summary_bsa.csv", row.names = F, quote=F)
write.csv(qtl_summary.bsa.round, "tables/qtl_summary_bsa_ROUNDED.csv", row.names = F, quote=F)
write.csv(qtl_summary.lm, "tables/qtl_summary_lm.csv", row.names = F, quote=F)
write.csv(qtl_summary.lm.round, "tables/qtl_summary_lm_ROUNDED.csv", row.names = F, quote=F)

