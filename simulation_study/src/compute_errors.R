args = commandArgs(trailingOnly=TRUE)

# The number of the dataset, one of 1-8 or 10
ds.idx <- as.integer(args[1])
# Number of MCMC iterations that were run, must match a value used in run_pseudo_MCMC.R or run_pseudo_MCMC_slurm.R
n.iter <- as.integer(args[2])

# In case treess is installed somewhere other than standard locations
if ( length(args) == 3 ) {
  .libPaths(args[3])
}

library(treess)

# Seed
source("simulation_study/src/utils.R")
set.seed(getSeed(dsnames[ds.idx],n.iter,base=47))

# Determine simulation sizes
n.samps <- 1000
n.thin <- n.iter/n.samps

DS <- list.files("simulation_study/output/",full.names=TRUE)

ds <- DS[grepl(paste0(dsnames[ds.idx],"_MCMC_ngen_",n.iter,"_samplefreq_",n.thin),DS,fixed=TRUE)]

# Get the pre-run MCMC, called "sims"
load(ds)

# Monte Carlo error, brute forced
brute.force.mce <- bruteForceMCMCSE(sims)

# ESS-based SE
rootedRF <- function(trees) {
  phangorn::RF.dist(trees,rooted=TRUE)
}

library(rrnni)
rNNIDist <- function(trees) {
  recover()
  ranked <- lapply(trees,as_ranked)
  
}
equivalent.mce <- effectiveSizeEquivalentError(sims,tree.dist=rootedRF,ess.methods=getESSMethods(recommended=TRUE))

out.file <- paste0(dsnames[ds.idx],"_mcse.Rdata")
save(brute.force.mce,equivalent.mce,file=out.file,version=2)

