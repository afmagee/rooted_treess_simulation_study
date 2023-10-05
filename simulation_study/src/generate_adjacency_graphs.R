library(treess)

DS <- list.files("simulation_study/data/",full.names=TRUE)
DS <- DS[grepl("best_trees",DS)]

for (ds in DS) {
  trprob <- readTreeProbs(ds,CI.width=1,type="simple")
  adjacency.graph <- constructAdjacencyGraph(trprob$trees,trprob$probs,rooted=TRUE,trim=TRUE)
  
  out.file <- gsub("data","output",ds)
  out.file <- gsub("_best_trees.txt","_adjacency_graph.Rdata",out.file,fixed=TRUE)
  
  save(adjacency.graph,file=out.file)
}
