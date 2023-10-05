library(treess)
library(coda)
library(parallel)

source("Normal_random_walk/src/simpleMVRW.R")

set.seed(42)

# Define sample sizes
ngen <- round(seq(1e3,1e5,length.out=200))
nsamp <- 1000

# How many replicates per sample size?
nrep <- 100

RMCE <- mclapply(ngen,function(n){
  thin <- n/nsamp
  x <- lapply(1:nrep,function(xx){
    simpleMVRW(n,0.3,1)[seq(thin,n,thin),]
  })
  ess <- unlist(lapply(x,effectiveSize))
  
  mu_hat <- mean(unlist(x))
  squared_deviations <- unlist(lapply(x,function(xx){
    (mean(xx) - mu_hat)^2
  }))
  mcmcse <- sqrt(mean(squared_deviations))
  
  equiv_x <- lapply(ess,function(neff) {
    rnorm(round(neff))
  })
  equiv_mu_hat <- mean(unlist(equiv_x))
  equiv_squared_deviations <- unlist(lapply(equiv_x,function(xx){
    (mean(xx) - equiv_mu_hat)^2
  }))
  essse <- sqrt(mean(equiv_squared_deviations))
  
  res <- c((mcmcse - essse)/mcmcse,mean(ess))
  names(res) <- c("RMCE","ESS")
  return(res)
},mc.cores=5,mc.preschedule=FALSE)
  
RMCE <- do.call(rbind,RMCE)
ESS.RMCE <- RMCE[,"ESS"]
RMCE <- RMCE[,"RMCE"]

save(RMCE,ESS.RMCE,file="Normal_random_walk/output/RMCE.Rdata",version=2)
