library(treess)
library(phangorn)
source("simulation_study/src/read_gzipped_trees.R")

DS <- list.files("simulation_study/data",full.names=TRUE)

for (ds in DS) {
  
  cat(ds,"\n")
  
  # Malagasy example datasets are stored as tar.gz files
  trees <- read.gz.nexus.trees(ds)
  
  # Count trees to get probabilities of topologies
  counted <- countTrees(trees,rooted=TRUE)
  
  # truncate to 4096 best trees
  probs <- counted$counts/sum(counted$counts)
  prob_cutoff <- min(which(cumsum(probs) >= 0.95)) # if there are none, this return Inf and the cutoff will be 4096 trees
  cutoff <- ifelse(prob_cutoff > 4096, 4096, prob_cutoff)
  if ( length(trees) < cutoff ) {
    cutoff <- length(trees)
  }
  
  trees <- counted$trees[1:cutoff]
  probs <- probs[1:cutoff]

  # write
  out.file <- gsub(".tar.gz","_best_trees.txt",ds)
  cat("",sep="",file=out.file,append=TRUE)
  for (i in 1:length(probs)) {
    pure_topo <- trees[[i]]
    pure_topo$edge.length <- NULL
    cat(probs[i]," ",write.tree(pure_topo),"\n",sep="",file=out.file,append=TRUE)
  }
}
