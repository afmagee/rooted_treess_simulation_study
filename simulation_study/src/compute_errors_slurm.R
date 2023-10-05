# In case treess is installed somewhere other than standard locations
new.lib.path <- ""
.libPaths(new.lib.path)

library(treess)
source("simulation_study/src/utils.R")

# Determine simulation sizes
n.iter <- 10^(3:7)
n.samps <- 1000

DS <- list.files("simulation_study/output/",full.names=TRUE)
DS <- DS[grepl("_MCMC_",DS)]

# slurm job ID
task.id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

slurm.n.iter <- rep(n.iter,length(dsnames))
slurm.dsnames <- rep(dsnames,length(n.iter))

start <- 3*(task.id-1)+1
end <- 3*task.id

for (i in start:end) {
  this.n.iter <- slurm.n.iter[i]
  this.n.thin <- this.n.iter/n.samps
  this.dsname <- slurm.dsnames[i]
  
  # Get the pre-run MCMC, called "sims"
  ds <- DS[grepl(paste0(this.dsname,"_MCMC_ngen_",this.n.iter,"_samplefreq_",this.n.thin),DS,fixed=TRUE)]
  load(ds)

  # for reproducibility
  set.seed(getSeed(this.dsname,this.n.iter,47))

  # Monte Carlo error, brute forced
  brute.force.mce <- bruteForceMCMCSE(sims)

  # ESS-based SE
  equivalent.mce <- effectiveSizeEquivalentError(sims,tree.dist="RF")
  
  out.file <- paste0("simulation_study/output/",this.dsname,"_error_ngen_",this.n.iter,"_samplefreq_",this.n.thin,".Rdata")
  save(brute.force.mce,equivalent.mce,file=out.file,version=2)
}
