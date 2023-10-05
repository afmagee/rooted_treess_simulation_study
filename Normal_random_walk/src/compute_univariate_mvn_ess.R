source("Normal_random_walk/src/simpleMVRW.R")

set.seed(42)

# Define sample sizes
ngen <- round(seq(1e3,1e5,length.out=20))
nsamp <- 1000

# How many replicates per sample size?
nrep <- 5

I10 <- lapply(ngen,function(n){
  thin <- n/nsamp
  ess <- sapply(1:nrep,function(xx){
    coda::effectiveSize(simpleMVRW(n,0.3,10,proposal="univariate")[seq(thin,n,thin),])
  })
  return(ess)
})

save(I10,file="Normal_random_walk/output/I10.Rdata")

avg.avg <- unlist(lapply(I10,function(x){mean(as.numeric(x))}))

plot(ngen,avg.avg,type="l")
abline(lm(avg.avg ~ ngen),col="red")

summary(lm(avg.avg ~ ngen))


plot(ngen,avg.avg,ylim=range(unlist(I10)),type="l",xlab="# samples",ylab="per-variable ESS")
for (i in 1:length(ngen)) {
  y <- as.numeric(I10[[i]])
  points(x=rep(ngen[i],length(y)),y=y)
}
