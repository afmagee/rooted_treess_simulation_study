set.seed(42)

library(treess)
library(phangorn)
library(methods)
source("~/git_repos/tree_convergence_code/real_data_examples/src/utils.R")
source("~/git_repos/tree_convergence_code/MCMC_validation/src/general_plots.R")

cop <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/cophyline.tar.gz")
gep <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/gephyromantis.tar.gz")
het <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/heterixalus.tar.gz")
par <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/paroedura.tar.gz")
phe <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/phelsuma.tar.gz")
uro <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/uroplatus.tar.gz")

malagasy <- list(cop,gep,het,par,phe,uro)
names(malagasy) <- c("Cophyline","Gephyromantis","Heterixalus","Paroedura","Phelsuma","Uroplatus")

for (i in 1:length(malagasy)) {
  herps <- malagasy[[i]]
  treeprobs <- countTrees(herps,rooted=TRUE)
  treeprobs$probs <- treeprobs$counts/sum(treeprobs$counts)
  adjacency.graph <- constructAdjacencyGraph(trees=treeprobs$trees,weights=treeprobs$probs,trim=TRUE,rooted=TRUE)
  
  mcmc <- simulatePhylogeneticMCMC(adjacency.graph,ngen=1e5,nchains=100)
  
  est.tree.probs <- lapply(1:length(adjacency.graph$trees),function(i){
    apply(mcmc$indices,2,function(x){
      sum(x == i)/length(x)
    })
  })
  
  out.file <- paste0("MCMC_validation/figures/",names(malagasy)[i],".pdf")
  
  pdf(out.file,width=5,height=3)
  par(mai=c(0.8,0.8,0.5,0.01))
  scatterPlotWithErrors(adjacency.graph$probs,est.tree.probs,bars="CI",xlab="true probability",ylab="estimated probability",pch=4,lwd=2,log="xy")
  dev.off()
  
}