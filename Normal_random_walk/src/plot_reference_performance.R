library(treess)
library(coda)

source("Normal_random_walk/src/simpleMVRW.R")

# Define sample sizes
ngen <- round(seq(1e3,5e4,length.out=50))
nsamp <- 1000

# How many replicates per sample size?
nrep <- 100

load("Normal_random_walk/output/RMCE.Rdata")

plot(ngen,RMCE)

pdf("Normal_random_walk/figures/RMCE_on_normals.pdf",width=5,height=2.5)
  par(mfrow=c(1,2),mai=c(0.8,0.5,0.05,0.05),omi=c(0.01,0.25,0.01,0.01))
  hist(RMCE,xlab="RMCE",main="",ylab="",freq=FALSE,border=NA,col="grey70")
  mtext("Density",2,2.5)
  hist(1/(1 - RMCE),xlab="1/(1 - RMCE)",main="",ylab="",freq=FALSE,border=NA,col="grey70")
dev.off()