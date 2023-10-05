# Tests rooted version of DS1-8 and DS10 as in bitbucket.org/afmagee/tree_convergence_code/
# Assumes adjacency graphs for those DS, as unrooted trees, are available

# set.seed(42)
# 
# library(treess)
# library(phangorn)
# library(methods)
# source("~/git_repos/tree_convergence_code/real_data_examples/src/utils.R")
# 
# source("~/git_repos/tree_convergence_code/MCMC_validation/src/general_plots.R")
# 
# DS <- list.files("~/git_repos/tree_convergence_code/simulation_study/data",full.names=TRUE)
# 
# for (ds in DS) {
#   treeprobs <- readTreeProbs(ds,1,"simple")
#   treeprobs$trees <- ape::root.multiPhylo(treeprobs$trees,1,resolve.root=TRUE)
#   adjacency.graph <- constructAdjacencyGraph(trees=treeprobs$trees,weights=treeprobs$probs,trim=TRUE,rooted=TRUE)
#   
#   mcmc <- simulatePhylogeneticMCMC(adjacency.graph,ngen=1e5,nchains=100)
#   
#   est.tree.probs <- lapply(1:length(adjacency.graph$trees),function(i){
#     apply(mcmc$indices,2,function(x){
#       sum(x == i)/length(x)
#     })
#   })
#   
#   out.file <- paste0("~/git_repos/tree_convergence_code/MCMC_validation/rooted_figures/",gsub("_adjacency_graph.Rdata","",basename(ds)),".pdf")
#   
#   pdf(out.file,width=5,height=3)
#   par(mai=c(0.8,0.8,0.5,0.01))
#   scatterPlotWithErrors(adjacency.graph$probs,est.tree.probs,bars="CI",xlab="true probability",ylab="estimated probability",pch=4,lwd=2,log="xy")
#   dev.off()
#   
# }
# 
# rm(list=ls())
