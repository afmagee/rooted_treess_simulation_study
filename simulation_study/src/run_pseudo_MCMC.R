args = commandArgs(trailingOnly=TRUE)

# The number of the dataset in the vector dsnames declared in utils
ds.idx <- as.integer(args[1])
# Number of MCMC iterations to run, in the paper this is a number in 10^(3:7)
n.iter <- as.integer(args[2])

# In case treess is installed somewhere other than standard locations
if ( length(args) == 3 ) {
  .libPaths(args[3])
}

library(treess)

# Seed
source("simulation_study/src/utils.R")
set.seed(getSeed(dsnames[ds.idx],n.iter,base=42))

# Determine simulation sizes
n.samps <- 1000
n.thin <- n.iter/n.samps

DS <- list.files("simulation_study/output/",full.names=TRUE)
DS <- DS[grepl("adjacency_graph",DS)]

ds <- DS[grepl(paste0(dsnames[ds.idx],"_adjacency_graph"),DS)]

# Get the pre-computed adjacency graph
load(ds)

sims <- simulatePhylogeneticMCMC(adjacency.graph,ngen=n.iter,nchains=100,thin=n.thin)
out.file <- gsub("_adjacency_graph",paste0("_MCMC_ngen_",n.iter,"_samplefreq_",n.thin),ds)
save(sims,file=out.file)

