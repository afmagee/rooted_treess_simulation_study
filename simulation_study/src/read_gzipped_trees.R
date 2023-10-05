# Reads nexus trees out of a gzipped file
# Currently, tip names are lost in the process
# Currently, it is assumed that all lines with trees include character "[&R]" in them
read.gz.nexus.trees <- function(file) {
  con <- gzfile(file)
  tmp <- scan(con,sep="\n",what=character())
  close(con)
  is_tree <- grepl("[&R]",tmp,fixed=TRUE)
  return(read.tree(text=tmp[is_tree]))
}

# Reads nexus trees out of a gzipped file
# Currently, tip names are lost in the process
# Currently, it is assumed that all lines with trees include character "[&R]" in them
get.logPosterior.gz.nexus.trees <- function(file) {
  con <- gzfile(file)
  tmp <- scan(con,sep="\n",what=character())
  close(con)
  is_tree <- grepl("[&R]",tmp,fixed=TRUE)
  # Only lines with trees in them
  trace <- tmp[is_tree]
  trace <- unlist(lapply(trace, function(x){
    x <- strsplit(x,"] = [&R]",fixed=TRUE)[[1]][1]
    x <- strsplit(x,"[&lnP=",fixed=TRUE)[[1]][2]
    x <- as.numeric(x)
    return(x)
  }))
  return(trace)
}

# tree STATE_20000 [&lnP=-60984.39420736647] = [&R]