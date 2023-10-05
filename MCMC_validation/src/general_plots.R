# plots measurements with errors, all of which have given true value
# vector x, list y
Ibar <- function(x,probs=c(0.05,0.95),true.val,mean.col="black",...) {
  # recover()
  
  ylim <- range(unlist(x))
  if ( hasArg("ylim") ) {
    ylim <- list(...)$ylim
  }
  bar_width <- 0.1
  if ( hasArg(bar.width) ) {
    bar_width <- list(...)$bar.width
  }
  the_order <- order(unlist(lapply(x,mean)),decreasing=TRUE)
  
  n <- length(x)
  plot(NULL,NULL,xlim=c(0.5,n+0.5),ylim=ylim,...)
  
  is_covered <- unlist(lapply(x,quantile,probs=min(probs))) < true.val & unlist(lapply(x,quantile,probs=max(probs))) > true.val
  
  bar_col <- rep("firebrick1",n)
  bar_col[is_covered] <- "springgreen4"
  bar_col <- adjustcolor(bar_col,0.5)
  
  # bar (vertical)
  for (i in 1:n) {
    lines(x=c(i,i),y=quantile(x[[the_order[i]]],probs=probs),col=bar_col[i],...)
  }
  
  # bar (cross)
  for (i in 1:n) {
    lines(x=c(i-bar_width/2,i+bar_width/2),y=quantile(x[[the_order[i]]],probs=probs[c(1,1)]),col=bar_col[i],...)
    lines(x=c(i-bar_width/2,i+bar_width/2),y=quantile(x[[the_order[i]]],probs=probs[c(2,2)]),col=bar_col[i],...)
  }
  
  # means
  points(unlist(lapply(x,mean))[the_order],...)
 
  abline(h=true.val,lty=2,col="grey",lwd=2)
  
}

# plots measurements with errors (y) against values in x
# vector x, list y
scatterPlotWithErrors <- function(x,y,bars=c("SE","CI"),probs=c(0.05,0.95),mean.col="black",one.to.one.col="grey",...) {
  # recover()
  
  n <- length(x)
  
  ymeans <- unlist(lapply(y,mean))
  
  if ( bars == "CI" ) {
    ymins <- unlist(lapply(y,quantile,prob=probs[1]))
    ymaxs <- unlist(lapply(y,quantile,prob=probs[2]))
  } else if ( bars == "SE" ) {
    warning("Should not use SE for MCMC samples, sqrt(ESS) is needed.")
    yr <- unlist(lapply(y,sd))
    yr <- yr/sqrt(lengths(y))
    ymins <- ymeans - yr/2
    ymaxs <- ymeans + yr/2
  } else {
    stop("Valid options for bars are SE|CI")
  }
  
  is_covered <- ymins < x & ymaxs > x
  
  bar_col <- rep("firebrick1",length(ymins))
  bar_col[is_covered] <- "springgreen4"
  bar_col <- adjustcolor(bar_col,0.5)
  
  ylim <- range(c(ymins,ymaxs))
  if ( hasArg("ylim") ) {
    ylim <- list(...)$ylim
  } else if ( hasArg("log") && ylim[1] == 0 ) {
    ylim[1] <- 0.1 * min(x)
    ymins[ymins == 0] <- ylim[1]
  }
  xlim <- ylim
  
  # start the plot
  plot(NULL,NULL,xlim=xlim,ylim=ylim,...)
  
  # deal with bar widths, especially if there's a log-scale issue
  bar_width <- 0.0075 * (par("usr")[2] - par("usr")[1])
  if ( hasArg(bar.width) ) {
    bar_width <- list(...)$bar.width
  }
  
  Ileft <- x - bar_width/2
  Iright <- x + bar_width/2
  
  if ( hasArg("log") ) {
    effective_widths <- (Iright - Ileft)/x
    equalize <- effective_widths[1]/effective_widths
    Ileft <- x - equalize * bar_width/2
    Iright <- x + equalize * bar_width/2
  }
  
  # 1:1 line
  abline(a=0,b=1,col=one.to.one.col,...)
  
  # bar (vertical)
  for (i in 1:n) {
    lines(x=c(x[i],x[i]),y=c(ymins[i],ymaxs[i]),col=bar_col,...)
  }
  
  # bar (cross)
  for (i in 1:n) {
    lines(x=c(Ileft[i],Iright[i]),y=c(ymins[i],ymins[i]),col=bar_col[i],...)
    lines(x=c(Ileft[i],Iright[i]),y=c(ymaxs[i],ymaxs[i]),col=bar_col[i],...)
  }
  
  # means
  points(x,ymeans,...)
  
}

slidingWindows <- function(x,y,n.windows) {
  breaks <- seq(min(x),max(x),length.out=n.windows+1)
  means <- sapply(1:n.windows,function(i){
    if (i == n.windows) {
      mean(y[x >= breaks[i] & x <= breaks[i+1]])
    } else {
      mean(y[x >= breaks[i] & x < breaks[i+1]])
    }
  })
  xw <- c(breaks[1],sort(rep(breaks[2:n.windows],2)),breaks[n.windows+1])
  yw <- means[sort(rep(1:n.windows,2))]
  return(list(x=xw,y=yw))
}

plotScaleBar <- function(x,heatcolors,n.labels) {
  plot(NULL,NULL,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
  txt <- round(min(x) + (max(x) - min(x)) * c(0:(n.labels-1))/(n.labels-1))
  txtwidth <- max(strwidth(txt))
  right <- (1 - txtwidth) * 0.975
  for (i in 1:length(heatcols)) {
    polygon(x=c(0,right,right,0),y=c(i-1,i-1,i,i)/length(heatcols),border=NA,col=heatcols[i])
  }
  y <- (seq(1,length(heatcols),length.out=n.labels) - 0.5)/length(heatcols)
  text(x=right+txtwidth/2,y=y,labels=txt)
}
