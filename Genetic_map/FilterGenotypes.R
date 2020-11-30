setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/genetic_map/")

##################################         Import and clean up data        #####################################

geno_file <- read.table("data/SquashF2_geno.txt", header = T)

#Make separate files for quality scores and genotypes
snps <- geno_file[,1:12]
geno <- geno_file[,13:ncol(geno_file)]

#Get rid of X in beginning of column names
colnames(geno) <- gsub(pattern="^X", replacement = "", x = colnames(geno))

##################################        Filter genotypes        #####################################

#Look at some stats
hist(snps$READ_DEPTH, breaks=50)
hist(snps$CALL_RATE)
hist(snps$MAF)

#First filter for just really high confidence markers
confidence_filter <- snps$MAF > 0.10 &
  snps$CALL_RATE > 90 &
  snps$READ_DEPTH > 15 &
  snps$READ_DEPTH < 90
geno <- geno[confidence_filter,]
snps <- snps[confidence_filter,]

#Look at individual call rate and get rid of high missing individuals
missing <- apply(geno,2, function(x) sum(is.na(x))/length(x))
hist(missing, breaks=100)
sum(missing>0.4)
colnames(geno)[missing>0.4]
geno <- geno[,missing<0.4]

#Turn into 0,1,2
geno<- geno*2

#Make separate data frames for parents and F2 individuals
geno_parent <- geno[,colnames(geno) %in% c("Dunja", "15_6015")]
geno_f2 <- geno[,!colnames(geno) %in% c("Dunja", "15_6015")]

#Now filter based on parental alleles and segregation

#Filter just snps that are fixed for opposite alleles in parents
parent_filter <- !is.na(geno_parent$Dunja) &
                 !is.na(geno_parent$'15_6015') &
                 ((geno_parent$Dunja == 0 & geno_parent$`15_6015` == 2) |
                    (geno_parent$Dunja == 2 & geno_parent$`15_6015` == 0))

geno_f2 <- geno_f2[parent_filter,]
geno_parent <- geno_parent[parent_filter,]
snps <- snps[parent_filter,]

#chisquared function for expected 1:2:1 segregation
chisq_test <- function(x){
  obs <- c(sum(x==0, na.rm=T), sum(x==1, na.rm=T), sum(x==2, na.rm=T))
  exp <- sum(obs)*c(0.25,.5,.25)
  chi2 <- sum((obs-exp)^2/exp)
  return(pchisq(chi2,2,lower.tail = F))
}
chisq <- apply(geno_f2, 1, chisq_test)
hist(chisq)
cols <- rep("red", nrow(geno_f2))
cols[as.integer(gsub("Cp4.1LG","",snps$CHROM))%%2==1] <- "blue"
plot(-log10(chisq), col=cols)
abline(h=2)
View(cbind(as.character(snps$CHROM), snps$BP, -log10(chisq), t(apply(geno_f2,1,table))))

#Filter based on chisq p-value < .01
chisq_filter <- chisq > .01
geno_f2 <- geno_f2[chisq_filter,]
geno_parent <- geno_parent[chisq_filter,]
snps <- snps[chisq_filter,]

#Remove markers on LG 00
chrom_filter <- snps$CHROM!="Cp4.1LG00"
snps <- snps[chrom_filter,]
snps$CHROM <- droplevels(snps$CHROM)
geno_f2 <- geno_f2[chrom_filter,]
geno_parent <- geno_parent[chrom_filter,]

############################### Find duplicated individuals ######################

#Identify possible duplicated individuals
concordance <- function(x,y){
  return(sum(x==y, na.rm=T)/sum(!is.na(x) & !is.na(y)))
}

d <- ncol(geno_f2)
pairwise_concordance <- matrix(NA, ncol=d, nrow=d)
for(i in 1:(d-1)){
  for(j in (i+1):d){
    pairwise_concordance[i,j] <- concordance(geno_f2[i], geno_f2[j])
  }
}

pdf("plots/concordance.pdf")
hist(pairwise_concordance, main = '', xlab = "% Pairwise Concordant Genotypes")
rug(pairwise_concordance)
dev.off()

clones <- which(pairwise_concordance > 0.9, arr.ind=T)
clone_names <- matrix(NA, nrow=nrow(clones), ncol=ncol(clones))
for(i in 1:nrow(clones)){
  for(j in 1:ncol(clones)){
    clone_names[i,j] <- colnames(geno_f2)[clones[i,j]]
  }
}

#There are 4 pairs of duplicated samples

#What is the average genotyping error rate among these duplicated samples?
1-mean(pairwise_concordance[clones])

#For each one we will remove the genotype with the least missing data
#And delete both of the phenotypic records later
clones_keep <- rep(NA, 4)
clones_remove <- rep(NA, 4)
for(i in 1:4){
  missing <- apply(geno_f2[,clones[i,]], 2, function(x) sum(is.na(x)))
  clones_keep[i] <- clone_names[i,which(missing == min(missing))]
  clones_remove[i] <- clone_names[i,which(missing == max(missing))]
}

geno_f2 <- geno_f2[,!colnames(geno_f2) %in% clones_remove]

############################### Remove redundant markers ######################

#Function to see if two markers have redundant information at all non-missing sites
dup <- function(x,y){
  return(all(x == y, na.rm=T))
}

#This function takes a group of redundant markers and get the ones to remove based on missingness
find_removals <- function(geno, coordinates){
  na_counts <- apply(geno[coordinates,], 1, function(x) sum(is.na(x)))
  best_markers <- coordinates[na_counts == min(na_counts)]
  if(length(best_markers) == 1){
    sampled_marker <- best_markers
  }else{
    sampled_marker <- sample(best_markers, 1)
  }
  return(coordinates[coordinates != sampled_marker])
}

# Go through markers and test duplicated to adjacent marker
# For groups of adjacent markers that are duplicated, just retain the one with the least missing data

i <- 1
markers_to_remove <- c()
d <- nrow(geno_f2)

set.seed(1957)
while(i < d){
  if(dup(geno_f2[i,], geno_f2[i+1,])){
    dup_run <- c(i,i+1)
    j <- i + 2
    run = T
    while(run){
      if(j > d){ #End when you get to the end of the markers
        bad_columns <- c(markers_to_remove, find_removals(geno_f2, dup_run))
        run = F
        i <- i + 2
      }else if(dup(geno_f2[i,], geno_f2[j,])){
        dup_run <- c(dup_run,j)
        j <- j + 1        
      }else{
        run = F
        i <- i + length(dup_run) - 1
        markers_to_remove <- c(markers_to_remove, find_removals(geno_f2, dup_run))
      }
    }
  }else{
    i <- i + 1
  }
}

#########################        Write filtered genotype and snp files ####################

geno_f2_prune <- geno_f2[-markers_to_remove,]
geno_parent_prune <- geno_parent[-markers_to_remove,]
snps_prune <- snps[-markers_to_remove,]

geno_filter_unite <- cbind(geno_f2, geno_parent)
geno_filter_unite_prune <- cbind(geno_f2_prune, geno_parent_prune)

# Number of markers per chromosome
table(snps$CHROM)
table(snps_prune$CHROM)

write.csv(geno_filter_unite, "tables/geno_filter.csv", quote=F, row.names = F)
write.csv(geno_filter_unite_prune, "tables/geno_filter_prune.csv", quote=F, row.names = F)
write.csv(snps, "tables/snps_filter.csv", quote=F, row.names=F)
write.csv(snps_prune, "tables/snps_filter_prune.csv", quote=F, row.names=F)

##############################        Make file for import into r/qtl        #####################################

#First turn 0,1,2 into ABH

A <- geno_parent_prune$Dunja
B <- geno_parent_prune$`15_6015`
for(i in 1:ncol(geno_f2_prune)){
  f2 <- geno_f2_prune[,i]
  f2[(f2 == 0 & A == 0) | (f2 == 2 & A == 2)] <- "A"
  f2[(f2 == 0 & B == 0) | (f2 == 2 & B == 2)] <- "B"
  geno_f2_prune[,i] <- f2
}
geno_f2_prune[geno_f2_prune == 1] <- "H"

abh_file <- cbind("CHROM" = snps_prune$CHROM, "BP" = snps_prune$BP, geno_f2_prune)
write.csv(abh_file, "tables/geno_prune_abh.csv", quote=F, row.names = F)

#Transpose genos
geno_f2_prune <- t(geno_f2_prune)

#Import phenotypes
blues <- read.csv("../F23/tables/f23_blues.csv")
colnames(blues)[1] <- "Sample" #Rename first column

blues$Sample[!blues$Sample %in% rownames(geno_f2_prune)]
#10 samples with only pheno data - including the 3 checks
rownames(geno_f2_prune)[!rownames(geno_f2_prune) %in% blues$Sample]
#1 sample with only geno data

#Get phenos in geno order
blues.f2 <- blues[match(rownames(geno_f2_prune), blues$Sample),]
#Replace the duplicated samples phenotypes with NA
blues[blues$Sample %in% clones_keep,2:6] <- NA
blues.f2[blues.f2$Sample %in% clones_keep,2:6] <- NA

#Get chromosomes as numbers
chrom_numbers <- as.integer(gsub("Cp4.1LG", "", snps_prune$CHROM))

#Make marker names
markers <- paste(chrom_numbers, snps_prune$BP, sep="_")

#Put it all together
header <- as.character(c(colnames(blues.f2)[2:ncol(blues.f2)], "Sample", markers))
r1 <- as.character(c(rep("", ncol(blues.f2)), chrom_numbers))
body <- as.matrix(cbind(blues.f2[,2:ncol(blues.f2)], rownames(geno_f2_prune), geno_f2_prune))
rqtl_geno <- rbind(r1,body)
colnames(rqtl_geno) <- header

#Save file
write.csv(rqtl_geno, "tables/f2_geno_rqtl.csv", quote=F, row.names=F)
#Save pheno file with NAs replacing the mismatched families
write.csv(blues, "tables/blues_final.csv", quote=F, row.names=F)

