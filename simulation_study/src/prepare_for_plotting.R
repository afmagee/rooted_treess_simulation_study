library(treess)
source("simulation_study/src/utils.R")

DS <- list.files("simulation_study/output/",full.names=TRUE)
DS <- DS[grepl("error",DS)]

# we ignore all tree/split split probs smaller than this
min.split.prob <- 0.01
min.tree.prob <- 1e-7

# loop over datasets, compute the MCMCSE and ESS-SE for all split/tree probabilities and MRC tree
monte.carlo.error <- lapply(DS,function(ds){
  # get Monte Carlo error
  load(ds)
  
  # Get DS and run info
  txt <- strsplit(ds,"/")[[1]]
  txt <- txt[length(txt)]
  txt <- gsub(".Rdata","",txt,fixed=TRUE)
  
  this.ds <- gsub("_","",strsplit(txt,"_",fixed=TRUE)[[1]][1],fixed=TRUE)
  this.ds <- as.numeric(gsub("DS","",this.ds))
  
  txt <- strsplit(txt,"ngen_")[[1]][2]
  ngen <- as.numeric(gsub("_","",strsplit(txt,"_",fixed=TRUE)[[1]][1],fixed=TRUE))
  
  # get "true" (aka best estimate) probs of trees/splits
  probs.file <- gsub("error","probs",ds)
  load(probs.file)
  
  # for getting the tree size
  mcmc.file <- gsub("error","MCMC",ds)
  load(mcmc.file)
  ntaxa <- length(sims$trees[[1]]$tip.label)
  
  # Brute-force "true" values
  per.split.var <- rowMeans(brute.force.mce$splitProbSquaredError)
  per.tree.var <- rowMeans(brute.force.mce$treeProbSquaredError)
  
  per.split.sd <- sqrt(per.split.var)
  per.tree.sd <- sqrt(per.tree.var)
  
  # truncate
  use.splits <- (all.probs$split.probs > min.split.prob) & (all.probs$split.probs < 1)
  use.trees <- (all.probs$tree.probs > min.tree.prob) & (all.probs$tree.probs < 1)
  
  true.per.split.sd <- per.split.sd[use.splits]
  true.per.tree.sd <- per.tree.sd[use.trees]
  
  true.split.probs <- all.probs$split.probs[use.splits]
  true.tree.probs <- all.probs$tree.probs[use.trees]

  true.ctsd <- sqrt(mean(brute.force.mce$MRCSquaredError))

  true.vals <- c(true.per.split.sd,true.per.tree.sd,true.ctsd)
  
  # loop over each ESS method to compute ESS-SE
  per.ess <- lapply(equivalent.mce,function(ess) {
    per.split.var <- rowMeans(ess$splitProbSquaredError)
    per.tree.var <- rowMeans(ess$treeProbSquaredError)

    per.split.sd <- sqrt(per.split.var)
    per.tree.sd <- sqrt(per.tree.var)
    
    equiv.ctsd <- sqrt(mean(ess$MRCSquaredError))
    
    # truncate
    per.split.sd <- per.split.sd[use.splits]
    per.tree.sd <- per.tree.sd[use.trees]
    
    res <- c(true.per.split.sd - per.split.sd,
             true.per.split.sd,
             true.split.probs,
             true.per.tree.sd - per.tree.sd,
             true.per.tree.sd,
             true.tree.probs,
             true.ctsd - equiv.ctsd,
             true.ctsd,
             mean(ess$ESS),
             sd(ess$ESS),
             this.ds,
             ngen)
    
    names(res) <- c(paste0("delta.split.sd.",1:length(per.split.sd)),
                    paste0("true.split.sd.",1:length(per.split.sd)),
                    paste0("split.prob.",1:length(per.split.sd)),
                    paste0("delta.tree.sd.",1:length(per.tree.sd)),
                    paste0("true.tree.sd.",1:length(per.tree.sd)),
                    paste0("tree.prob.",1:length(per.tree.sd)),
                    "delta.ctsd",
                    "true.ctsd",
                    "ESS",
                    "SD.ESS",
                    "DS",
                    "ngen")
    return(res)
  })
  
  return(per.ess)
})

save(monte.carlo.error,file="simulation_study/output/MCMC_and_ESS_SE.Rdata")

rm(list=ls())