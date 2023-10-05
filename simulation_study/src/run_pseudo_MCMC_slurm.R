# In case treess is installed somewhere other than standard locations
new.lib.path <- ""
.libPaths(new.lib.path)

library(treess)
source("simulation_study/src/utils.R")

# Determine simulation sizes
n.iter <- 10^(3:7)
n.samps <- 1000
n.thin <- n.iter/n.samps

DS <- list.files("simulation_study/output/",full.names=TRUE)
DS <- DS[grepl("adjacency_graph",DS)]

# slurm job ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Get the pre-computed adjacency graph
this.ds <- dsnames[touse]
ds <- DS[grepl(this.ds,DS)]
load(ds)

for (i in 1:length(n.iter)) {
  # for reproducibility
  set.seed(getSeed(this.ds,n.iter,42))
  
  # run and store fake MCMC
  sims <- simulatePhylogeneticMCMC(adjacency.graph,ngen=n.iter[i],nchains=100,thin=n.thin[i])
  out.file <- gsub("_adjacency_graph",paste0("_MCMC_ngen_",n.iter[i],"_samplefreq_",n.thin[i]),ds)
  save(sims,file=out.file,version=2)
}
