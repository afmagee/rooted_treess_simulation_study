library(latex2exp)
library(treess)
source("simulation_study/src/utils.R")

# reference performance
load("Normal_random_walk/output/RMCE.Rdata")

ess.cutoffs <- c(250,500)

############
# Loop over ESS threshold for trees and splits
############
for (ess.thresh in ess.cutoffs ) {
  # Since we overwrite a few names a few times, we need to re-do a few things in each loop
  # \begin{obnoxiousReloading}
  
  # names of all methods
  all.methods <- c(getESSMethods(),"logPosteriorESS")
  
  # Processed results
  if (exists("monte.carlo.error")) {
    rm(monte.carlo.error)
  }
  load("simulation_study/output/MCMC_and_ESS_SE.Rdata")
  # \end{obnoxiousReloading}
  
  # construct shell of list for split prob plot
  split.probs.by.ess <- list(low=list(RMCE=vector("list",length(all.methods)),
                                      ITMCE=vector("list",length(all.methods)),
                                      probs=vector("list",length(all.methods))),
                             high=list(RMCE=vector("list",length(all.methods)),
                                       ITMCE=vector("list",length(all.methods)),
                                       probs=vector("list",length(all.methods)))
                         )
  
  names(split.probs.by.ess$low$RMCE) <- all.methods
  names(split.probs.by.ess$low$ITMCE) <- all.methods
  names(split.probs.by.ess$low$probs) <- all.methods
  names(split.probs.by.ess$high$RMCE) <- all.methods
  names(split.probs.by.ess$high$ITMCE) <- all.methods
  names(split.probs.by.ess$high$probs) <- all.methods
  
  for (i in 1:45) {
    for (j in 1:12) {
      this.ess <- monte.carlo.error[[i]][[j]]["ESS"]
      this.method <- names(monte.carlo.error[[i]])[j]
      # For fixedN, we split into ngen <= 1e4 and ngen > 1e4 instead of splitting on ESS
      if ( this.method == "fixedN" ) {
        this.ess <- ifelse(monte.carlo.error[[i]][[j]]["ngen"] <= 1e4,ess.thresh-1,ess.thresh+1)
      }
      metrics <- names(monte.carlo.error[[i]][[j]])
      
      # Compute RMCE and ITMCE for all splits for this ESS in this simulation
      these.rmce <- monte.carlo.error[[i]][[j]][grepl("delta.split.sd",metrics)]/monte.carlo.error[[i]][[j]][grepl("true.split.sd",metrics)]
      these.itmce <- 1 / (1 - these.rmce)
      these.probs <- monte.carlo.error[[i]][[j]][grepl("split.prob",metrics)]
      
      # Assign RMCE and ITMCE to low or high performance regimes
      if (this.ess < ess.thresh) {
        split.probs.by.ess$low$RMCE[[this.method]] <- c(split.probs.by.ess$low$RMCE[[this.method]],these.rmce)
        split.probs.by.ess$low$ITMCE[[this.method]] <- c(split.probs.by.ess$low$ITMCE[[this.method]],these.itmce)
        split.probs.by.ess$low$probs[[this.method]] <- c(split.probs.by.ess$low$probs[[this.method]],these.probs)
      } else {
        split.probs.by.ess$high$RMCE[[this.method]] <- c(split.probs.by.ess$high$RMCE[[this.method]],these.rmce)
        split.probs.by.ess$high$ITMCE[[this.method]] <- c(split.probs.by.ess$high$ITMCE[[this.method]],these.itmce)
        split.probs.by.ess$high$probs[[this.method]] <- c(split.probs.by.ess$high$probs[[this.method]],these.probs)
      }
    }
  }
  
  # construct shell of list for tree prob plot
  tree.probs.by.ess <- list(low=list(RMCE=vector("list",length(all.methods)),
                                      ITMCE=vector("list",length(all.methods)),
                                      probs=vector("list",length(all.methods))),
                             high=list(RMCE=vector("list",length(all.methods)),
                                       ITMCE=vector("list",length(all.methods)),
                                       probs=vector("list",length(all.methods)))
  )
  
  names(tree.probs.by.ess$low$RMCE) <- all.methods
  names(tree.probs.by.ess$low$ITMCE) <- all.methods
  names(tree.probs.by.ess$low$probs) <- all.methods
  names(tree.probs.by.ess$high$RMCE) <- all.methods
  names(tree.probs.by.ess$high$ITMCE) <- all.methods
  names(tree.probs.by.ess$high$probs) <- all.methods
  
  for (i in 1:45) {
    for (j in 1:12) {
      this.ess <- monte.carlo.error[[i]][[j]]["ESS"]
      this.method <- names(monte.carlo.error[[i]])[j]
      # For fixedN, we tree into ngen <= 1e4 and ngen > 1e4 instead of treeting on ESS
      if ( this.method == "fixedN" ) {
        this.ess <- ifelse(monte.carlo.error[[i]][[j]]["ngen"] <= 1e4,ess.thresh-1,ess.thresh+1)
      }
      metrics <- names(monte.carlo.error[[i]][[j]])
      
      # Compute RMCE and ITMCE for all trees for this ESS in this simulation
      these.rmce <- monte.carlo.error[[i]][[j]][grepl("delta.tree.sd",metrics)]/monte.carlo.error[[i]][[j]][grepl("true.tree.sd",metrics)]
      these.itmce <- 1 / (1 - these.rmce)
      these.probs <- monte.carlo.error[[i]][[j]][grepl("tree.prob",metrics)]
      
      # Assign RMCE and ITMCE to low or high performance regimes
      if (this.ess < ess.thresh) {
        tree.probs.by.ess$low$RMCE[[this.method]] <- c(tree.probs.by.ess$low$RMCE[[this.method]],these.rmce)
        tree.probs.by.ess$low$ITMCE[[this.method]] <- c(tree.probs.by.ess$low$ITMCE[[this.method]],these.itmce)
        tree.probs.by.ess$low$probs[[this.method]] <- c(tree.probs.by.ess$low$probs[[this.method]],these.probs)
      } else {
        tree.probs.by.ess$high$RMCE[[this.method]] <- c(tree.probs.by.ess$high$RMCE[[this.method]],these.rmce)
        tree.probs.by.ess$high$ITMCE[[this.method]] <- c(tree.probs.by.ess$high$ITMCE[[this.method]],these.itmce)
        tree.probs.by.ess$high$probs[[this.method]] <- c(tree.probs.by.ess$high$probs[[this.method]],these.probs)
      }
    }
  }
  
  # # construct shell of list for MRC tree plot
  # ctsd.by.ess <- list(low=list(RMCE=vector("list",length(all.methods)),
  #                              ITMCE=vector("list",length(all.methods))),
  #                     high=list(RMCE=vector("list",length(all.methods)),
  #                               ITMCE=vector("list",length(all.methods)))
  #                     )
  # 
  # names(ctsd.by.ess$low$RMCE) <- all.methods
  # names(ctsd.by.ess$low$ITMCE) <- all.methods
  # names(ctsd.by.ess$high$RMCE) <- all.methods
  # names(ctsd.by.ess$high$ITMCE) <- all.methods
  # 
  # for (i in 1:45) {
  #   for (j in 1:12) {
  #     this.ess <- monte.carlo.error[[i]][[j]]["ESS"]
  #     this.method <- names(monte.carlo.error[[i]])[j]
  #     # For fixedN, we split into ngen <= 1e4 and ngen > 1e4 instead of splitting on ESS
  #     if ( this.method == "fixedN" ) {
  #       this.ess <- ifelse(monte.carlo.error[[i]][[j]]["ngen"] <= 1e4,ess.thresh-1,ess.thresh+1)
  #     }
  #     metrics <- names(monte.carlo.error[[i]][[j]])
  #     
  #     # Compute RMCE and ITMCE
  #     # For RMCE, if the true SD is 0 and we estimate the SD to be 0, we define the RMCE to be 0
  #     these.rmce <- monte.carlo.error[[i]][[j]][grepl("delta.ctsd",metrics)]/monte.carlo.error[[i]][[j]][grepl("true.ctsd",metrics)]
  #     these.rmce[monte.carlo.error[[i]][[j]][grepl("delta.ctsd",metrics)] == 0 & monte.carlo.error[[i]][[j]][grepl("true.ctsd",metrics)] == 0] <- 0
  #     these.itmce <- 1 / (1 - these.rmce)
  #     
  #     if (this.ess < ess.thresh) {
  #       ctsd.by.ess$low$RMCE[[this.method]] <- c(ctsd.by.ess$low$RMCE[[this.method]],these.rmce)
  #       ctsd.by.ess$low$ITMCE[[this.method]] <- c(ctsd.by.ess$low$ITMCE[[this.method]],these.itmce)
  #     } else {
  #       ctsd.by.ess$high$RMCE[[this.method]] <- c(ctsd.by.ess$high$RMCE[[this.method]],these.rmce)
  #       ctsd.by.ess$high$ITMCE[[this.method]] <- c(ctsd.by.ess$high$ITMCE[[this.method]],these.itmce)
  #     }
  #   }
  # }
  
  # Vertical ordering of plotted ESS methods
  plot.orders <- list(
    rev(c(5,8,9,1,3,12)),
    rev(c(5,10,8,9,4,11,2,1,6,7,3,12))
  )
  
  main.or.supp <- c("_main","_supplemental")
  
  metrics <- c("_RMCE","_ITMCE")
  
  heights <- c(2.25,6)
  
  # xlabs <- c(TeX("$\\frac{(SE_{MCMC} - SE_{ESS})}{SE_{MCMC}}$"),TeX("$\\frac{1}{1 - RMCE}$"))
  xlabs <- c("RMCE","ITMCE")
  
  for (edx in 1:2) {
    
    XLAB <- xlabs[edx]
    
    to.plot.lines <- c(quantile(RMCE,c(0.1,0.25)),0,quantile(RMCE,c(0.75,0.9)))
    if (edx == 2) {
      to.plot.lines <- c(quantile(1/(1-RMCE),c(0.1,0.25)),1,quantile(1/(1-RMCE),c(0.75,0.9)))
    }
    
    for (idx in 1:2) {
      
      toplot <- plot.orders[[idx]]
      
      OMI <- c(0.4,1.15,0.01,0.01)
      if ( idx == 2 ) {
        OMI <- c(0.4,2.1,0.01,0.01)
      }
      
      # cat("ess.thresh = ",ess.thresh,"; edx = ",edx,"; idx = ",idx,"\n",sep="")
      # cat(class(split.probs.by.ess$low[[edx]][["fullApproximateESS"]]),"\n")
      # cat(class(split.probs.by.ess$high[[edx]][["fullApproximateESS"]]),"\n")
      # cat("names(split.probs.by.ess) yields: ",paste0(names(split.probs.by.ess),collapse="; "),"\n",sep="")
      # cat("names(split.probs.by.ess$low) yields: ",paste0(names(split.probs.by.ess$low),collapse="; "),"\n",sep="")
      # cat("names(split.probs.by.ess$low[[edx]]) yields: ",paste0(names(split.probs.by.ess$low[[edx]]),collapse="; "),"\n",sep="")
      # cat("names(split.probs.by.ess$low[[edx]][toplot]) yields: ",paste0(names(split.probs.by.ess$low[[edx]][toplot]),collapse="; "),"\n",sep="")
      
      #########
      # splits
      #########
      
      # Some RMCE can be 1, so some ITMCE can be infinite (but we don't really care, that's an edge case where proportions are a bit silly)
      split.measures.finite <- unlist(split.probs.by.ess$low[[edx]][toplot])
      split.measures.finite <- split.measures.finite[is.finite(split.measures.finite)]
      
      best.split.measures.finite <- unlist(split.probs.by.ess$high[[edx]][c("minPseudoESS","frechetCorrelationESS")])
      best.split.measures.finite <- best.split.measures.finite[is.finite(best.split.measures.finite)]
      
      pdf(paste0("simulation_study/figures/split_performance",main.or.supp[idx],metrics[edx],"_ESS_cutoff_",ess.thresh,".pdf"),width=7,height=heights[idx])
        par(mfrow=c(1,3),mai=c(0.05,0.05,0.15,0.05),omi=OMI)
        plotPerformance(split.probs.by.ess$low[[edx]][toplot],
                        split.probs.by.ess$low$probs[toplot],
                        xlim=quantile(split.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"))
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=2,at=1:length(toplot),line=0.5,text=names(split.probs.by.ess$low[[edx]])[toplot],las=1,cex=0.7)
        mtext(side=3,text=paste0("ESS < ",ess.thresh),line=0)
        
        legend.coords <- c(-6,6)
        if (idx == 2) {
          legend.coords <- c(-30,12)
        }
        if (edx == 2) {
          legend.coords[1] <- "blah"
        }
        plotPerformance(split.probs.by.ess$high[[edx]][toplot],
                        split.probs.by.ess$high$probs[toplot],
                        xlim=quantile(split.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"),
                        legend.xy=legend.coords)
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
        
        plotPerformance(split.probs.by.ess$high[[edx]][toplot],
                        split.probs.by.ess$high$probs[toplot],
                        xlim=quantile(best.split.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"))
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
      dev.off()

      #########
      # tree probs
      #########
      
      # Some RMCE can be 1, so some ITMCE can be infinite (but we don't really care, that's an edge case where proportions are a bit silly)
      tree.measures.finite <- unlist(tree.probs.by.ess$low[[edx]][toplot])
      tree.measures.finite <- tree.measures.finite[is.finite(tree.measures.finite)]
      
      best.tree.measures.finite <- unlist(tree.probs.by.ess$high[[edx]][c("minPseudoESS","frechetCorrelationESS")])
      best.tree.measures.finite <- best.tree.measures.finite[is.finite(best.tree.measures.finite)]
      
      
      pdf(paste0("simulation_study/figures/tree_performance",main.or.supp[idx],metrics[edx],"_ESS_cutoff_",ess.thresh,".pdf"),width=7,height=heights[idx])
        par(mfrow=c(1,3),mai=c(0.05,0.05,0.15,0.05),omi=OMI)
        plotPerformance(tree.probs.by.ess$low[[edx]][toplot],
                        tree.probs.by.ess$low$probs[toplot],
                        scale.color.to.max.prob=TRUE,
                        log.scale.color=TRUE,
                        xlim=quantile(tree.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"),
                        max.points.per.row=1000)
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=2,at=1:length(toplot),line=0.5,text=names(tree.probs.by.ess$low[[edx]][toplot]),las=1,cex=0.7)
        mtext(side=3,text=paste0("ESS < ",ess.thresh),line=0)
        
        legend.coords <- c(-27.5,6)
        if (idx == 2) {
          legend.coords <- c(-35,12)
        }
        if (edx == 2) {
          legend.coords[1] <- "blah"
        }
        plotPerformance(tree.probs.by.ess$high[[edx]][toplot],
                        tree.probs.by.ess$high$probs[toplot],
                        scale.color.to.max.prob=TRUE,
                        log.scale.color=TRUE,
                        xlim=quantile(tree.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"),
                        max.points.per.row=1000,
                        legend.xy=legend.coords)
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
        
        plotPerformance(tree.probs.by.ess$high[[edx]][toplot],
                        tree.probs.by.ess$high$probs[toplot],
                        scale.color.to.max.prob=TRUE,
                        log.scale.color=TRUE,
                        xlim=quantile(best.tree.measures.finite,c(0.005,0.995)),
                        line.at=to.plot.lines,
                        line.lty=c(2,1,1,1,2),
                        line.lwd=c(1.25,1.25,2,1.25,1.25),
                        line.col=c("red","red","black","red","red"),
                        max.points.per.row=1000)
        mtext(side=1,line=2.25,text=XLAB,cex=0.9)
        mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
      dev.off()
      
      #########
      # summary tree
      #########
      
      # # Some RMCE can be -Inf for this, where the true SD is 0 but we estimate non-zero
      # mrc.measures.finite <- unlist(ctsd.by.ess$low[[edx]][toplot])
      # mrc.measures.finite <- mrc.measures.finite[is.finite(mrc.measures.finite)]
      # 
      # best.mrc.measures.finite <- unlist(ctsd.by.ess$high[[edx]][c("minPseudoESS","frechetCorrelationESS")])
      # best.mrc.measures.finite <- best.mrc.measures.finite[is.finite(best.mrc.measures.finite)]
      # 
      # pdf(paste0("simulation_study/figures/MRC_performance",main.or.supp[idx],metrics[edx],"_ESS_cutoff_",ess.thresh,".pdf"),width=7,height=heights[idx]*0.75)
      #   par(mfrow=c(1,3),mai=c(0.05,0.05,0.25,0.05),omi=OMI)
      #   plotPerformance(ctsd.by.ess$low[[edx]][toplot],
      #                   NA,
      #                   xlim=quantile(mrc.measures.finite,c(0.005,0.995)),
      #                   line.at=to.plot.lines,
      #                   line.lty=c(2,1,1,1,2),
      #                   line.lwd=c(1.25,1.25,2,1.25,1.25),
      #                   line.col=c("red","red","black","red","red"))
      #   mtext(side=1,line=2.25,text=XLAB,cex=0.9)
      #   mtext(side=2,at=1:length(toplot),line=0.5,text=names(ctsd.by.ess$low[[edx]][toplot]),las=1,cex=0.7)
      #   mtext(side=3,text=paste0("ESS < ",ess.thresh),line=0)
      #   
      #   plotPerformance(ctsd.by.ess$high[[edx]][toplot],
      #                   NA,
      #                   xlim=quantile(mrc.measures.finite,c(0.005,0.995)),
      #                   line.at=to.plot.lines,
      #                   line.lty=c(2,1,1,1,2),
      #                   line.lwd=c(1.25,1.25,2,1.25,1.25),
      #                   line.col=c("red","red","black","red","red"))
      #   mtext(side=1,line=2.25,text=XLAB,cex=0.9)
      #   mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
      #   
      #   plotPerformance(ctsd.by.ess$high[[edx]][toplot],
      #                   NA,
      #                   xlim=quantile(best.mrc.measures.finite,c(0.005,0.995)),
      #                   line.at=to.plot.lines,
      #                   line.lty=c(2,1,1,1,2),
      #                   line.lwd=c(1.25,1.25,2,1.25,1.25),
      #                   line.col=c("red","red","black","red","red"))
      #   mtext(side=1,line=2.25,text=XLAB,cex=0.9)
      #   mtext(side=3,text=paste0("ESS >= ",ess.thresh),line=0)
      # dev.off()
      
    }
  }
  
}

############
# Pool all analyses for summary trees
############

# Since we overwrite a few names a few times, we need to re-do a few things in each loop
# \begin{obnoxiousReloading}

# names of all methods
all.methods <- c(getESSMethods(),"logPosteriorESS")

# Processed results
if (exists("monte.carlo.error")) {
  rm(monte.carlo.error)
}
load("simulation_study/output/MCMC_and_ESS_SE.Rdata")
# \end{obnoxiousReloading}


# construct shell of list for MRC tree plot
ctsd.by.ess <- list(all=list(RMCE=vector("list",length(all.methods)),
                             ITMCE=vector("list",length(all.methods)),
                             ESS=vector("list",length(all.methods)))
)

names(ctsd.by.ess$all$RMCE) <- all.methods
names(ctsd.by.ess$all$ITMCE) <- all.methods
names(ctsd.by.ess$all$ESS) <- all.methods

for (i in 1:45) {
  for (j in 1:12) {
    this.ess <- monte.carlo.error[[i]][[j]]["ESS"]
    this.method <- names(monte.carlo.error[[i]])[j]
    metrics <- names(monte.carlo.error[[i]][[j]])
    
    # Compute RMCE and ITMCE
    # For RMCE, if the true SD is 0 and we estimate the SD to be 0, we define the RMCE to be 0
    these.rmce <- monte.carlo.error[[i]][[j]][grepl("delta.ctsd",metrics)]/monte.carlo.error[[i]][[j]][grepl("true.ctsd",metrics)]
    these.rmce[monte.carlo.error[[i]][[j]][grepl("delta.ctsd",metrics)] == 0 & monte.carlo.error[[i]][[j]][grepl("true.ctsd",metrics)] == 0] <- 0
    these.itmce <- 1 / (1 - these.rmce)
    
    ctsd.by.ess$all$RMCE[[this.method]] <- c(ctsd.by.ess$all$RMCE[[this.method]],these.rmce)
    ctsd.by.ess$all$ITMCE[[this.method]] <- c(ctsd.by.ess$all$ITMCE[[this.method]],these.itmce)
    ctsd.by.ess$all$ESS[[this.method]] <- c(ctsd.by.ess$all$ESS[[this.method]],this.ess)
  }
}

# Vertical ordering of plotted ESS methods
plot.orders <- list(
  rev(c(5,8,9,1,3,12)),
  rev(c(5,10,8,9,4,11,2,1,6,7,3,12))
)

main.or.supp <- c("_main","_supplemental")

metrics <- c("_RMCE","_ITMCE")

heights <- c(2.25,6)

# xlabs <- c(TeX("$\\frac{(SE_{MCMC} - SE_{ESS})}{SE_{MCMC}}$"),TeX("$\\frac{1}{1 - RMCE}$"))
xlabs <- c("RMCE","ITMCE")

for (edx in 1:2) {
  
  
  XLAB <- xlabs[edx]
  
  to.plot.lines <- c(quantile(RMCE,c(0.1,0.25)),0,quantile(RMCE,c(0.75,0.9)))
  if (edx == 2) {
    to.plot.lines <- c(quantile(1/(1-RMCE),c(0.1,0.25)),1,quantile(1/(1-RMCE),c(0.75,0.9)))
  }
  
  for (idx in 1:2) {
    
    legend.coords[1] <- "blah"
    if (idx == 1 && edx == 2) {
      legend.coords <- c(8,6)
    }
    if (idx == 2 && edx == 1) {
      legend.coords <- c(-8,12)
    }
    
    toplot <- plot.orders[[idx]]
    
    OMI <- c(0.6,1.2,0.01,0.01)
    if ( idx == 2 ) {
      OMI <- c(0.6,2.15,0.01,0.01)
    }
    
    # Some RMCE can be -Inf for this, where the true SD is 0 but we estimate non-zero
    mrc.measures.finite <- unlist(ctsd.by.ess$all[[edx]][toplot])
    mrc.measures.finite <- mrc.measures.finite[is.finite(mrc.measures.finite)]
    
    best.mrc.measures.finite <- unlist(ctsd.by.ess$high[[edx]][c("minPseudoESS","frechetCorrelationESS")])
    best.mrc.measures.finite <- best.mrc.measures.finite[is.finite(best.mrc.measures.finite)]
    
    pdf(paste0("simulation_study/figures/MRC_performance",main.or.supp[idx],metrics[edx],"_all_ESS.pdf"),width=5.5,height=heights[idx]*0.75)
      
      par(mfrow=c(1,2),mai=c(0.05,0.05,0.05,0.05),omi=OMI)
      plotPerformance(ctsd.by.ess$all[[edx]][toplot],
                      ctsd.by.ess$all$ESS[toplot],
                      xlim=quantile(mrc.measures.finite,c(0.005,0.995)),
                      line.at=to.plot.lines,
                      line.lty=c(2,1,1,1,2),
                      line.lwd=c(1.25,1.25,2,1.25,1.25),
                      line.col=c("red","red","black","red","red"),
                      # quantile.colors=c("#FFFFFF00","grey60"),
                      color.scheme="ESS",
                      legend.xy=legend.coords,
                      log.scale.color=TRUE)
      mtext(side=1,line=2.25,text=XLAB,cex=0.9)
      mtext(side=2,at=1:length(toplot),line=0.5,text=names(ctsd.by.ess$all[[edx]][toplot]),las=1,cex=0.7)
      
      xl <- c(-0.8,0.8)
      if ( edx == 2 ) {
        xl <- c(0.4,2.0)
      }
      
      plotPerformance(ctsd.by.ess$all[[edx]][toplot],
                      ctsd.by.ess$all$ESS[toplot],
                      xlim=xl,
                      line.at=to.plot.lines,
                      line.lty=c(2,1,1,1,2),
                      line.lwd=c(1.25,1.25,2,1.25,1.25),
                      line.col=c("red","red","black","red","red"),
                      # quantile.colors=c("#FFFFFF00","grey60"),
                      color.scheme="ESS",
                      legend.xy=legend.coords,
                      log.scale.color=TRUE)
      mtext(side=1,line=2.25,text=XLAB,cex=0.9)

    dev.off()
  }

}
