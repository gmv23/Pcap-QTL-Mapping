library(qtl)
library(parallel)

squash <- readRDS("squash_map.rds")

run_permutation <- function(i){
  perm.out <- scantwo(squash, pheno.col=5, n.perm=20,
                      method="hk", incl.markers=T)
  saveRDS(perm.out, file = paste("out/perms", i, ".rds", sep=""))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1918)
batches <- 1:50
mclapply(batches, run_permutation, mc.set.seed=T)

