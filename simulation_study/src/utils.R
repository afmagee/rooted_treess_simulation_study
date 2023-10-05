# For consistent retrieval of datasets by index
dsnames <- c("cophyline","gephyromantis","heterixalus","paroedura","phelsuma","uroplatus")

# For consistent seeds regardless of method of calling
getSeed <- function(dsname,n.iter,base=42) {
  # Safe for n.iter > 1000 and base < 100
  return(base+which(dsnames == dsname)+n.iter)
}

plotPerformance <- function(performance,
                            color.by,
                            scale.color.to.max.prob=FALSE,
                            log.scale.color=FALSE,
                            xlim=NA,
                            line.at=0,
                            line.lwd=1,
                            line.lty=1,
                            line.col="black",
                            plot.quantiles=TRUE,
                            quantile.colors=c("grey80","grey60"),
                            median.color="grey30",
                            max.points.per.row=NA,
                            legend.xy=NA,
                            color.scheme="probs") {
  
  # recover()
  color.scheme <- tolower(color.scheme)
  
  if ( is.na(xlim[1]) ) {
    xlim <- range(unlist(performance),na.rm=TRUE)
  }
  
  plot(NULL,NULL,ylim=c(0.5,length(performance)+0.5),xlim=xlim,xlab="",ylab="",yaxt="n")
  
  # For coloring points, if desired
  bmin <- NA
  bmax <- NA
  n <- 101
  breaks <- NA
  
  if ( color.scheme != "probs" ) {
    bmin <- min(unlist(color.by))
    bmax <- max(unlist(color.by))
    n <- 101
    breaks <- seq(bmin,bmax,length.out=n)
    if ( log.scale.color ) {
      bmin <- min(unlist(color.by))
      breaks <- exp(seq(log(bmin),log(bmax),length.out=n))
    }
  } else {
    bmin <- 0
    bmax <- ifelse(scale.color.to.max.prob,max(unlist(color.by)),1)
    n <- 101
    breaks <- seq(bmin,bmax,length.out=n)
    if ( log.scale.color ) {
      bmin <- min(unlist(color.by))
      breaks <- exp(seq(log(bmin),log(bmax),length.out=n))
    }
  }
  
  for (i in 1:length(performance)) {
    if ( !is.null(performance[[i]]) ) {
      
      x <- performance[[i]]
      y <- runif(length(x),i-0.25,i+0.25)
      
      if ( plot.quantiles ) {
        mx <- median(x)
        polygon(x=quantile(x,probs=c(0.1,0.9))[c(1,2,2,1)],y=i+c(-0.4,-0.4,+0.4,+0.4),border=NA,col=quantile.colors[1])
        polygon(x=quantile(x,probs=c(0.25,0.75))[c(1,2,2,1)],y=i+c(-0.475,-0.475,+0.475,+0.475),border=NA,col=quantile.colors[2])
      }
      
      idx <- 1:length(x)
      if ( is.numeric(max.points.per.row) && length(x) > max.points.per.row ) {
        # Avoid plotting all points, but bias towards including the higher-probability points
        idx1 <- order(color.by[[i]],decreasing=TRUE)[1:floor(max.points.per.row/3)]
        avail <- !(1:length(x) %in% idx1)
        idx2 <- sample((1:length(x))[avail],max.points.per.row - length(idx1),replace=FALSE,prob=color.by[[i]][avail])
        idx <- c(idx1,idx2)
      }
      
      cols <- rep("black",length(x))
      viridis_pal <- viridis::viridis(101)
      if ( color.scheme == "ess" ) {
        viridis_pal <- viridis::viridis(101,option="inferno")
      }
      if ( is.list(color.by) ) {
        cols <- cut(color.by[[i]],breaks,labels=FALSE)
        cols <- viridis_pal[cols]
        cols <- cols[idx]
        
        p <- color.by[[i]][idx]
        x <- x[idx]
        y <- y[idx]
        
        y <- sort(y)
        x <- x[order(p)]
        cols <- cols[order(p)]
        points(x,y,col=cols,pch=16,cex=0.6)
      } else {
        points(x[idx],y[idx],col=cols[idx],pch=16,cex=0.6)
      }
      
      if ( plot.quantiles ) {
        lines(x=rep(mx,2),y=c(i-0.35,i+0.35),col=median.color,lwd=3)
      }
    }  
  }
  
  if ( is.numeric(line.at) ) {
    abline(v=line.at,lwd=line.lwd,lty=line.lty,col=line.col)
  }
  
  if ( is.numeric(legend.xy) && length(legend.xy) == 2 ) {
    xspan <- abs(par("usr")[2] - par("usr")[1])
    yspan <- abs(legend.xy[2] - par("usr")[3])
    xc <- c(legend.xy[1],legend.xy[1]+0.05*xspan)
    ymax <- legend.xy[2]
    ymin <- ymax - 0.5*yspan
    yc <- seq(ymin,ymax,length.out=length(viridis_pal)+1)
    for ( i in 1:length(viridis_pal) ) {
      polygon(x=xc[c(1,2,2,1)],y=yc[c(i,i,i+1,i+1)],col=viridis_pal[i],border=FALSE)
    }
    # labs <- round(pretty(breaks),2)
    # labs_at <- sapply(labs,function(x){which.min( (breaks - x)^2 )})
    labs_at <- round(seq(1,101,length.out=5))
    cex_labs <- 1.0
    labs <- breaks[labs_at]
    if ( color.scheme == "ess" ) {
      labs <- round(labs)
      cex_labs <- 0.75
    } else {
      if ( log.scale.color ) {
        labs <- formatC(labs,format="e",digits=1)
      } else {
        labs <- round(labs,2)
      }
    }
    text(x=xc[2],y=yc[labs_at],labels=labs,adj=0,cex=cex_labs)
  }
  
}

# A ridge plot, but histograms >> KDEs
historidge <- function(x,
                       colors,
                       global.summary.line=NULL,
                       global.summary.line.color="#66666690",
                       global.summary.line.lty=2,
                       global.summary.line.lwd=1.5,
                       xlab=NULL,
                       ylab="",
                       xline=2.5,
                       yline=2,
                       names=NULL,
                       hist.summary=NULL,
                       lines.lwd=1,
                       lines.lty=1,
                       summary.line.color=NULL,
                       nbreaks=21,
                       cex.xlab=1,
                       cex.ylab=1,
                       cex.xaxis=1,
                       cex.yaxis=1,
                       plot.y.axis=TRUE,
                       plot.x.axis=TRUE,
                       box=TRUE,
                       main="",
                       xmin=NULL,
                       xmax=NULL,
                       scale=0.95,
                       xaxis.pad.percent=0.01) {
  # recover()
  
  if ( class(x) != "list" ) {
    stop("Input must be list.")
  }
  
  all_datapoints <- unlist(x)
  
  if (any(is.na(all_datapoints))) {
    warning("There are NAs or NaNs in data. Removing and continuing")
  }
  
  summary_range <- range(all_datapoints,na.rm=TRUE)
  
  # pad a little
  xlim <- summary_range
  xlim[1] <- xlim[1] * ifelse(xlim[1] < 0,1+xaxis.pad.percent,1/(1+xaxis.pad.percent))
  xlim[2] <- xlim[2] * ifelse(xlim[2] < 0,1/(1+xaxis.pad.percent),1+xaxis.pad.percent)
  
  if (is.numeric(xmin)) {
    xlim[1] <- xmin
  }
  
  if (is.numeric(xmax)) {
    xlim[2] <- xmax
  }
  
  breaks <- seq(summary_range[1],summary_range[2],length.out=nbreaks)
  
  # Calculate the histograms
  x_hist <- lapply(x, hist, plot=FALSE, breaks=breaks)
  
  # We want all histograms to be on the same y-axis scale
  # So we need to make sure that the bars are standardized to the max bar height
  x_max <- max(unlist(lapply(x_hist,function(xh){max(xh$density)})))
  for (i in 1:length(x)) {
    x_hist[[i]]$density <- x_hist[[i]]$density / x_max * scale
  }
  
  # Graphical parameters
  y_max <- length(x) - 1 + max(x_hist[[length(x)]]$density)
  line_max <- 1.125 # for plotting summary lines
  
  # We need to be able to plot outside the margins
  current_par <- par(no.readonly=TRUE)
  reset_xpd <- FALSE
  if ( !(current_par$xpd == TRUE || is.na(current_par$xpd)) ) {
    stop("Please set par(xpd=TRUE) or par(xpd=NA) to par(xpd=TRUE).")
    # stop("Overriding par(xpd=TRUE) or par(xpd=NA) to par(xpd=TRUE). This will be reset but may cause issues.")
    # old_xpd <- current_par$xpd
    # current_par$xpd <- TRUE
    # reset_xpd <- TRUE
    # par(current_par)
  }
  
  # Start new plot
  frame()
  plot.window(xlim=xlim,ylim=c(0,y_max*1.01),main=main)
  
  if ( box ) {
    box()
  }
  
  if ( plot.x.axis ) {
    axis(side=1,cex.axis=cex.xaxis)
  }
  
  if ( !is.null(xlab) && xlab != "" ) {
    mtext(text=xlab,side=1,line=xline,cex=cex.xlab)
  }
  
  
  if ( plot.y.axis ) {
    
    # if (scale > 1) {
    #   warning("Cannot plot densities for each sub-axis, plotting only 0s")
    #   axis(side=2,cex.axis=cex.yaxis,labels=rep(0,length(x)),at=1:length(x)-1)
    # } else {
    #   subaxis <- vector("list",length(x))
    #   ticks_at <- vector("list",length(x))
    #   for (i in 1:length(x)) {
    #     subaxis[[i]] <- seq(0,max(x_hist[[i]]$density),length.out=3)
    #     ticks_at[[i]] <- i - 1 + seq(0,scale,length.out=3)
    #   }
    #   subaxis <- round(unlist(subaxis),1)
    #   ticks_at <- unlist(ticks_at)
    #   axis(side=2,cex.axis=cex.yaxis,labels=subaxis,at=ticks_at)
    # }
    
    if ( !is.null(names) && names != "" ) {
      # Label histograms
      box_width <- par("usr")[2] - par("usr")[1]
      box_height <- par("usr")[4] - par("usr")[3]
      y_tick_length <- box_width/box_height * strheight("hellodave") * par("tcl")
      tick_left = par("usr")[1] + y_tick_length
      # text(x=rep(tick_left,length(x)),y=1:length(x)-1,labels=names,cex=cex.xaxis,pos=2)
      text(x=rep(par("usr")[1],length(x)),y=1:length(x)-1,labels=names,cex=cex.xaxis,pos=2)
    }
    
  }
  
  if ( !is.na(ylab) && ylab != "") {
    mtext(text=ylab,side=2,line=yline,cex=cex.ylab)
  }
  
  # histograms
  for (i in length(x):1) {
    for (j in 1:(nbreaks-1)) {
      polygon(y=i - 1 + c(0,0,x_hist[[i]]$density[j],x_hist[[i]]$density[j]),
              x=c(x_hist[[i]]$breaks[j],x_hist[[i]]$breaks[j+1],x_hist[[i]]$breaks[j+1],x_hist[[i]]$breaks[j]),
              border=NA,col=colors[i])
    }
    if (!is.null(hist.summary)) {
      x <- hist.summary(x[[i]])
      for (j in 1:length(x)) {
        lines(y=c(i - line_max,i),x=rep(y[j],2),col=summary.line.color,lwd=lines.lwd,lty=lines.lty)
      }
    }
  }
  
  if ( is.numeric(global.summary.line) ) {
    lines(x=rep(global.summary.line,2),y=c(0,y_max*1.01),col=global.summary.line.color,lty=global.summary.line.lty,lwd=global.summary.line.lwd)
  }
  
  # if ( reset_xpd ) {
  #   current_par$xpd <- old_xpd
  #   par(current_par)
  # }
  # 
}