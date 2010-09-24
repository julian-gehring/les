######################################################################
## exported functions for les package
######################################################################

######################################################################
## Les constructor
######################################################################
Les <- function(pos, pval, chr)  {

  ## default values
  if(missing(chr))
    chr <- rep(0, length(pos))

  ## check inputs
  if(length(pos) != length(pval)  || length(pos) != length(chr))
    stop("'pos', 'pval' and 'chr' must have the same length.")
  if(any(is.na(pos)))
    stop("'pos' must not contain NAs.")
  if(any(pos %% 1 != 0))
    stop("'pos' must be a vector of integers.")
  if(!is.numeric(pval) || min(pval, na.rm=TRUE) < 0 || max(pval, na.rm=TRUE) > 1)
    stop("'pval' must be in the range [0,1].")

  ## throw out NAs in pval
  indValid <- !is.na(pval)
  pos <- as.integer(pos[indValid])
  chr <- factor(chr[indValid])

  ## sort
  ord <- order(chr, pos)
  pos <- pos[ord]
  pval <- pval[ord]
  chr <- chr[ord]

  object <- new(Class="Les",
                pos=pos, pval=pval, chr=chr, nChr=nlevels(chr),
                state="Les")
  
  return(object)
}


######################################################################
## triangWeight
######################################################################
triangWeight <- function(distance, win)  {
  
  weight <- 1 - abs(distance)/win
    
  return(weight)
}

######################################################################
## gaussWeight
######################################################################
gaussWeight <- function(distance, win)  {

  weight <- stats::dnorm(distance, sd=win/2)

  return(weight)
}


######################################################################
## rectangWeight
######################################################################
rectangWeight <- function(distance, win)  {

  n <- length(distance)
  weight <- rep(1/n, n)

  return(weight)
}


######################################################################
## epWeight
######################################################################
epWeight <- function(distance, win)  {
  
  weight <- 0.75*(1 - (distance/win)^2)
    
  return(weight)
}



######################################################################
## quartWeight
######################################################################
quartWeight <- function(distance, win)  {
  
  weight <- 15/16*(1 - (distance/win)^2)^2
    
  return(weight)
}


######################################################################
## tricubeWeight
######################################################################
tricubeWeight <- function(distance, win)  {
  
  weight <- 35/32*(1 - (distance/win)^2)^3
    
  return(weight)
}
