##################################################
## calcSingle
##################################################
calcSingle <- function(ind0, pos, pval, win,
                       weighting, grenander, se,
                       nBoot, conf, minProbes, custom)  {

  ## cut down
  pos0 <- pos[ind0]
  indCut <- seq(ind0-win, ind0+win)
  posCut <- pos[indCut]
  pvalCut <- pval[indCut]
  
  distance <- posCut - pos0
  indValid <- abs(distance) <= win
  # nValidProbes <- length(unique(pvalCut[indValid]))
  nValidProbes <- length(pvalCut[indValid])
  nUniqueProbes <- length(unique(pvalCut[indValid]))

  if(nUniqueProbes > minProbes)  {
    ## if enough probes
    dis <- distance[indValid]
    weight <- weighting(dis, win)
    indWeight <- weight > 0
    nWeight <- sum(indWeight)

    ## apply
    if(nBoot == FALSE)  {
      res <- fitGsri(pvalCut[indValid][indWeight], index=NULL,
                     weight[indWeight], nWeight, grenander, se,
                     custom)
    }
    else  {
      bo <- boot::boot(pvalCut[indValid][indWeight], fitGsri, nBoot,
                 cweight=weight[indWeight], nValidProbes=nWeight,
                 grenander=grenander, se=se, custom=custom)
      ## try not needed? TO DO !
      suppressWarnings(ci <- try(boot::boot.ci(bo, conf, type="perc"),
                                 silent=TRUE))
      if(class(ci) == "try-error")  ## TO DO: why this?
        res <- c(NA, NA)
      else
        res <- ci$perc[4:5]
    }
  }
  else  {
    ## if not enough probes
    if(nBoot == FALSE)
      res <- c(NA, NA, nValidProbes)
    else
      res <- c(NA, NA)
  }
  
  return(res)
}
## ok ## nValidProbes as input?


##################################################
## fitGsri
##################################################
fitGsri <- function(pval, index=NULL, cweight,
                    nValidProbes, grenander, se, custom)  {
  
  noBoot <- is.null(index)
  if(noBoot == TRUE)
    cdf <- wcdf2(pval, cweight, grenander)
  else
    cdf <- wcdf2(pval[index], cweight[index], grenander)
  if(any(cdf$cdf < 0))
    stop("weights < 0")
  res <- itLinReg(cdf$pval, cdf$cdf, cweight, nValidProbes, se, custom, noBoot)
  
  return(res)
}


##################################################
## itLinReg
##################################################
itLinReg <- function(x, y, cweight, nValidProbes, se, custom, noBoot)  {
  browser()
  maxIter <- nValidProbes
  restOld <- 0
  rest <- restOld
  q <- 1
  x <- x - 1
  y <- y - 1
  if(custom == TRUE)
    cw <- diagSquare(cweight, nValidProbes)
  ## iterative fitting
  for(i in 1:maxIter)  {
    rest <- nValidProbes - ceiling(q*nValidProbes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(c(nValidProbes-1, rest))
    if(is.na(rest) || restOld == rest)
      break
    ind <- rest:nValidProbes
    if(custom == TRUE)  {
      q <- slopeWeight(x[ind], y[ind], cw[ind,ind])
    }
    else  {
      xi <- x[ind]
      dim(xi) <- c(length(ind), 1)
      q <- qrSlope(xi, y[ind], cweight[ind])
    }
    restOld <- rest
  }

  ## return values
  if(noBoot == TRUE)  {
    if(se == TRUE && length(unique(x)) > 1)
      ses <- seFast(x, y, q)
    else
      ses <- NA
    res <- c(1 - min(q, 1), ses, nValidProbes)
  }
  else  {
    res <- 1 - min(q, 1, na.rm=TRUE)
  }
  
  return(res)
}


##################################################
## wcdf2
##################################################
wcdf2 <- function(pval, weight, grenander=FALSE)  {
  
  ord <- sort.list(pval, method="quick", na.last=NA)
  pvalSort <- pval[ord]
  weightSort <- weight[ord]
  nPval <- length(pvalSort)

  if(any(duplicated(pvalSort)))  {
    tabCount <- rle(pvalSort)$lengths  ## table()
    tabWeight <- table(pvalSort, weightSort)
    weightSum <- tabWeight %*% sort(unique(weightSort))
    cdf1 <- cumsum(weightSum)
    cdf1r <- rep.int(cdf1, tabCount)
    cdf2r <- rep.int(c(0, cdf1[-length(cdf1)]), tabCount)
    cdf2r <- cdf2r/cdf1r[nPval]
    cdf1r <- cdf1r/cdf1r[nPval]
    cdf <- cdf1r - (cdf1r - cdf2r)/2
    # weight2 <- rep.int(cdf, tabCount)
  }
  else  {
    cdf <- cumsum(weightSort)
    cdf <- cdf/cdf[nPval]
    cdf <- cdf - (cdf-c(0, cdf[-nPval]))/2
  }

  if(grenander == TRUE)
    cdf <- GSRI:::grenanderInterp(pvalSort, cdf)
  
  res <- list(pval=pvalSort, cdf=cdf)

  return(res)
  }


##################################################
## seFast
##################################################
seFast <- function(x, y, b)  {
  
  si     <- sum((y-b*x)^2)/(length(x)-1)
  se     <- as.numeric(sqrt(si/(t(x)%*%x)))
  
  return(se)
}
## ok ##




##################################################
## mcapply
##################################################
mcsapply <- function(X, FUN, ..., mc.cores=NULL)  {

  mcLoaded <- any(.packages() %in% "multicore") &&
  any("mclapply" %in% objects("package:multicore"))
  
  if(mcLoaded == TRUE && !is.null(mc.cores))  {
    res <- multicore::mclapply(X, FUN, ..., mc.cores=mc.cores)
    res <- sapply(res, c)
  }
  else  {
    res <- sapply(X, FUN, ...)
  }

  return(res)
}
## ok ##



##################################################
## gsri
##################################################
gsri <- function(pval, grenander=FALSE, se=TRUE, custom=FALSE)  {

  cweight <- rep(1, length(pval))
  res <- fitGsri(pval, NULL, cweight, length(pval),
                 grenander=grenander, se=se, custom=custom)
  res <- c(res, res[1]*res[3])
  names(res) <- c("GSRI", "se", "n", "nReg")

  return(res)
}


log2ind <- function(log)  {

  ind <- seq(along=log)[log]
  ind <- as.integer(ind)

  return(ind)
}
## ok ##


ind2log <- function(ind, n)  {

  log <- logical(n)
  log[ind] <- TRUE

  return(log)
}
## ok ##



##################################################
## slopeWeight
##################################################
slopeWeight <- function(x, y, c)  {
  s <- t(x) %*% c  ## do outside and subset? possible since matrix? no !
  b <- ((s %*% y) / (s %*% x))[1]  ## same as as.numeric()
  return(b)
}


qrSlope <- function(x, y, w)  {
  b <- as.numeric(lm.wfit(x, y, w)$coefficients)
  return(b)
}


diagSquare <- function(w, n)  {
  y <- rep(0, n*n)
  y[1+0:(n-1)*(n+1)] <- w
  dim(y) <- c(n, n)
  return(y)
}

