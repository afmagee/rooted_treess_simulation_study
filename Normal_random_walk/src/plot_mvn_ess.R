# Define sample sizes
ngen <- round(seq(1e3,1e5,length.out=20))
nsamp <- 1000

# How many replicates per sample size?
nrep <- 5

load("Normal_random_walk/output/M10.Rdata")

all.ess <- do.call(rbind,lapply(M10,function(x){do.call(rbind,x)}))

per.ess <- sapply(colnames(all.ess),function(m){
  sapply(1:length(ngen),function(i) {
    first <- (i-1)*nrep+1
    last <- i*nrep
    mean(all.ess[first:last,m])
  })
})

pdf("Normal_random_walk/figures/mvESS_vs_n.pdf",width=6,height=4)
  par(mai=c(1,1,0.01,0.01))
  plot(NULL,NULL,xlim=c(1e3,1e5),ylim=range(as.numeric(per.ess)),log="",xlab="chain length",ylab="multivariate ESS")
  for (i in 1:dim(per.ess)[2]) {
    lines(ngen,per.ess[,i],col=ifelse(i == 10,"red","black"))
  }
dev.off()
