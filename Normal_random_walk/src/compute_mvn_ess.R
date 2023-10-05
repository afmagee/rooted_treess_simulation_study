library(treess)

source("Normal_random_walk/src/simpleMVRW.R")

set.seed(42)

# distance-based treESS measures we can use on MVNs, ignore "fixedN" too
applicable.methods <- getESSMethods()[-c(3,10)]

# Define sample sizes
ngen <- round(seq(1e3,1e5,length.out=20))
nsamp <- 1000

# How many replicates per sample size?
nrep <- 5

M10 <- lapply(ngen,function(n){
  thin <- n/nsamp
  x <- lapply(1:nrep,function(xx){
    simpleMVRW(n,0.3,10)[seq(thin,n,thin),]
  })
  ess <- treess(x,dist,methods=applicable.methods)
  # Add multivariate ESS of MCMCSE as comparison
  for (i in 1:nrep) {
    ess[[i]]["multiESS"] <- mcmcse::multiESS(x[[i]])
  }
  return(ess)
})

if ( !dir.exists("Normal_random_walk/output") ) {
  dir.create("Normal_random_walk/output")
}

save(M10,file="Normal_random_walk/output/M10.Rdata")
