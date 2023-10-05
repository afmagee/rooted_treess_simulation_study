# Targets a diagonal multivariate normal with all-0 mean-vector and all variances 1
# Can either jointly propose moves to all dimensions, or loop over each dimension separately (in which case the chain is actually run for ngen*dimension steps and sampled every dimension steps)
simpleMVRW <- function(ngen,prop.sd,dimension,proposal="joint") {
  mcmc <- matrix(nrow=ngen,ncol=dimension)
  x <- rnorm(dimension,0,1)
  for (i in 1:ngen) {
    if ( proposal == "joint" ) {
      x_prop <- x + rnorm(dimension,0,prop.sd)
      log_AR <- sum(dnorm(x_prop,0,1,log=TRUE)) - sum(dnorm(x,0,1,log=TRUE))
      if ( log(runif(1)) < log_AR ) {
        x <- x_prop
      }
    } else if ( proposal == "univariate" ) {
      for (j in 1:dimension) {
        x_prop <- x
        x_prop[j] <- x_prop[j] + rnorm(1,0,prop.sd)
        log_AR <- sum(dnorm(x_prop,0,1,log=TRUE)) - sum(dnorm(x,0,1,log=TRUE))
        if ( log(runif(1)) < log_AR ) {
          x <- x_prop
        }
      }
    }
    mcmc[i,] <- x
  }
  return(mcmc)
}
