set.seed(42)

library(treess)
library(phangorn)

source("MCMC_validation/src/general_plots.R")

## Uniform probabilities

# Here our target is all 4-tip rooted trees
# All trees have equal probability
fake.ds.adj <- constructAdjacencyGraph(allTrees(4,TRUE),rooted=TRUE)

mcmc <- simulatePhylogeneticMCMC(fake.ds.adj,ngen=1e4,nchains=100)

# summarize and plot
est.tree.probs <- lapply(1:length(fake.ds.adj$trees),function(i){
  apply(mcmc$indices,2,function(x){
    sum(x == i)/length(x)
  })
})

pdf("MCMC_validation/figures/4_taxon_uniform.pdf",width=5,height=3)
  par(mai=c(0.75,0.8,0.5,0.01))
  Ibar(est.tree.probs,true.val=1/length(fake.ds.adj$probs),pch=4,main="4 taxon uniform distribution",xlab="",xaxt="n",ylab="estimated probability",lwd=2)
  axis(1,at=1:length(est.tree.probs),labels=1:length(est.tree.probs))
  mtext("Tree index",1,line=2)
dev.off()



## Non-uniform probabilities

# Here our target is all 5-tip rooted trees
# Probabilities are dirichlet-distributed
fake.ds.adj$probs <- sort(rgamma(length(fake.ds.adj$trees),1,1),decreasing=TRUE)
fake.ds.adj$probs <- fake.ds.adj$probs/sum(fake.ds.adj$probs)

mcmc <- simulatePhylogeneticMCMC(fake.ds.adj,ngen=1e4,nchains=100)

# summarize and plot
est.tree.probs <- lapply(1:length(fake.ds.adj$trees),function(i){
  apply(mcmc$indices,2,function(x){
    sum(x == i)/length(x)
  })
})

pdf("MCMC_validation/figures/4_taxon_dirichlet.pdf",width=5,height=3)
  par(mai=c(0.8,0.8,0.5,0.01))
  scatterPlotWithErrors(fake.ds.adj$probs,est.tree.probs,bars="CI",xlab="true probability",ylab="estimated probability",pch=4,lwd=2)
dev.off()

