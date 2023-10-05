# # In case treess is installed somewhere other than standard locations
# .libPaths("path/to/treess")

library(treess)
library(parallel)
source("simulation_study/src/utils.R")


DS <- list.files("simulation_study/output/",full.names=TRUE)
DS <- DS[grepl("MCMC",DS)]

# loop over datasets
monte.carlo.probs <- parallel::mclapply(DS,function(ds){
  
  # Get simulated trees
  load(ds)
  
  # Get tree coordinates
  coords <- as.RFcoords(sims$trees)$coords
  
  ntrees <- length(sims$trees)
  all.idx <- as.integer(sims$indices)
  tree.probs <- sapply(1:length(sims$trees),function(k){
    sum(all.idx == k)/length(all.idx)
  })
  
  split.probs <- colSums(coords * tree.probs)
  
  all.probs <- list(split.probs=split.probs,tree.probs=tree.probs)
  
  out.file <- gsub("MCMC","probs",ds)
  
  save(all.probs,file=out.file,version=2)
  
},mc.cores=15)
  