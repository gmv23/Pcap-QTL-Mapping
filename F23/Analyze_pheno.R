setwd("~/Documents/work/Smart_lab/P_capsici/QTL_mapping/F23/")
library(lme4)
library(RColorBrewer)
##################################         Import and clean up data        #####################################

disease <- read.csv("data/f23_ratings.csv")
disease$Rep <- as.factor(disease$Rep)

#Get rid of phenotypes that we don't have on all individuals
disease <- disease[,!colnames(disease) %in% c("dpi4", "dpi6", "dpi12")]

#Look at distribution of Number plants germed per plot
hist(disease$Germ)
sum(disease$Germ < 6)

#Get rid of plots that have fewer than 6 plants
disease <- disease[disease$Germ > 5,]

#Get rid of extra check "18_G_105"
disease <- disease[disease$entry!="18_G_105",]
disease$entry <- droplevels(disease$entry)

######################################    Get H2, BLUEs, checks, medians    #######################################

#Also make plot with distribution of residuals

#Make separate data frame without checks
checks <- c("17_G_020", "Dunja", "17_G_14")
disease_prog <- disease[!disease$entry %in% checks,]
disease_prog$entry <- droplevels(disease_prog$entry)
trait_columns <- 6:10
results <- data.frame("Trait" = colnames(disease)[trait_columns],
                      "H2" = NA,
                      "Dunja" = NA,
                      "Pc-NY21" = NA,
                      "F1" = NA,
                      "F23_RAN" = NA,
                      "F23_SEL" = NA)

trait_blues <- matrix(NA, nrow = nlevels(disease$entry), ncol = 5)
rownames(trait_blues) <- levels(disease$entry)
colnames(trait_blues) <- colnames(disease)[6:10]

#Which were 'random' F23s vs 'selected' F23s
selected_families <- grep("18_G", rownames(trait_blues))
random_families <- grep("Vog", rownames(trait_blues))

#Get harmonic mean of number of reps (very slightly unbalanced data)
r <- table(disease_prog$entry)
n_harm <- length(r)/sum(1/r)

#Make plot of residuals distributions for each trait
pdf("plots/residual_distributions.pdf", height = 10, width = 6)
old.par <- par(no.readonly = T)
par(mfrow=c(5,2))

#Loop through traits
for(i in 1:5){
  trait_name <- colnames(trait_blues)[i]

  #Fit random effects model
  mod1 <- lmer(disease_prog[,trait_columns[i]] ~ (1|disease_prog$entry) + (1|disease_prog$Rep))
  
  #Get H2
  trait.var <- as.data.frame(VarCorr(mod1))
  vg <- trait.var$vcov[1]
  ve <- trait.var$vcov[3]
  h2 <- vg/(vg + ve/n_harm)

  #Fit fixed effects model
  pheno <- disease[,trait_columns[i]]
  mod2 <- lmer(pheno ~ 1 + disease$entry + (1|disease$Rep))
 
  #Plot histogram of residuals as well as residuals vs fitted values
  hist(residuals(mod2), main=trait_name, xlab = "Residuals")
  plot(fitted(mod2), residuals(mod2), main=trait_name, xlab="Fitted values", ylab="Residuals")
  
  #Get BLUEs - add intercept to all entries, make intercept the first entry
  effects <- fixef(mod2)
  blues <- effects + effects[1]
  blues[1] <- effects[1]
  names(blues) <- levels(disease$entry)
  trait_blues[,i] <- blues[match(names(blues), rownames(trait_blues))]
  
  #Pull out check means, plus medians of selected vs non-selected families
  dunja <- blues["Dunja"]
  pcny21 <- blues["17_G_14"]
  f1 <- blues["17_G_020"]
  ran_median <- median(blues[random_families])
  sel_median <- median(blues[selected_families])
  
  #Add to results
  results[i,2:ncol(results)] <- round(c(h2, dunja, pcny21, f1, ran_median, sel_median),2)
}
par(old.par)
dev.off()

trait_blues <- as.data.frame(trait_blues)

#Save data
write.csv(results, "tables/trait_summaries.csv", quote=F)
write.csv(trait_blues, "tables/f23_blues.csv", quote=F)
######################################    Show phenotypic distribution    #######################################

pdf("plots/audpc_hist.pdf", width=4, height = 3)
old.par <- par(no.readonly = T)
par(mar=c(5,4,2,1), oma=c(1,1,1,1))
m <- rbind(c(1,1,1,2))
layout(m)
hist(trait_blues[c(random_families, selected_families),5],
     xlim = c(0,100),
     ylim=c(0,40),
     breaks=20,
     col=colors()[431],
     xlab = "rAUDPC BLUEs",
     ylab = "Frequency",
     main = "")
abline(v=dunja, lwd=2.5, lty = 1)
abline(v=f1, lwd=2.5, lty=2)
abline(v=pcny21, lwd=2.5, lty=3)
plot(0,type='n', xaxt='n', yaxt='n', ylab='', xlab='', bty='n', xlim=c(0,10), ylim=c(0,10))
legend(-20,9,
       legend = c(expression("Dunja F"[1]), "PcNY-21", expression("F"[1])),
       lty=c(1,3,2), bty='n', cex=1, lwd=2, xpd=NA,
       y.intersp=2)
dev.off()



