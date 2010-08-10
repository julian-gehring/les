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
      res <- les:::fitGsri(pvalCut[indValid][indWeight], index=NULL,
                           weight[indWeight], nWeight, grenander, se,
                           custom)
    }
    else  {
      bo <- boot::boot(pvalCut[indValid][indWeight], les:::fitGsri, nBoot,
                       cweight=weight[indWeight], nValidProbes=nWeight,
                       grenander=grenander, se=se, custom=custom)
      suppressWarnings(ci <- try(boot::boot.ci(bo, conf, type="perc"),
                                 silent=TRUE))
      if(class(ci) == "try-error")
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


##################################################
## fitGsri
##################################################
fitGsri <- function(pval, index=NULL, cweight,
                    nValidProbes, grenander, se, custom)  {

  noBoot <- is.null(index)
  if(noBoot == TRUE)  {
    cdf <- les:::wcdf2(pval, cweight, FALSE)
  }
  else  {
    cdf <- les:::wcdf2(pval[index], cweight[index], FALSE)
  }
  res <- les:::itLinReg(cdf$pval, cdf$cdf, cweight, nValidProbes, se, custom, noBoot)

  if(grenander == TRUE)  {
    cdf$cdf <- les:::cdfCorrect(cdf$pval, cdf$cdf, 1-res[1])
    cdf$cdf <- les:::grenanderPass(cdf$pval, cdf$cdf, cdf$unique)
    cdf$cdf <- GSRI:::grenanderInterp(cdf$pval, cdf$cdf)
    res <- les:::itLinReg(cdf$pval, cdf$cdf, cweight, nValidProbes, se, custom, noBoot)
  }
  
  return(res)
}


##################################################
## cdfCorrect
##################################################
cdfCorrect <- function(x, y, q0)  {

  z <- y
  indLower <- y < q0*x
  indUpper <- y > 1 - q0*(1 - x)
  z[indLower] <- q0*x[indLower]
  z[indUpper] <- 1 - q0*(1 - x[indUpper])

  return(z)
}


##################################################
## itLinReg
##################################################
itLinReg <- function(x, y, cweight, nValidProbes, se, custom, noBoot)  {

  maxIter <- nValidProbes
  rEstOld <- 0
  rEst <- rEstOld
  q <- 1
  x <- x - 1
  y <- y - 1
  if(custom == TRUE)  {
    cw <- les:::diagSquare(cweight, nValidProbes)
  }
  ## iterative fitting
  for(i in 1:maxIter)  {
    rEst <- nValidProbes - ceiling(q*nValidProbes)
    rEst <- max(c(rEstOld, rEst, 1))
    rEst <- min(c(nValidProbes, rEst))
    if(is.na(rEst) || rEstOld == rEst)
      break
    ind <- rEst:nValidProbes
    if(custom == TRUE)  {
      q <- les:::slopeWeight(x[ind], y[ind], cw[ind,ind])
    }
    else  {
      xi <- x[ind]
      dim(xi) <- c(length(ind), 1)
      q <- les:::qrSlope(xi, y[ind], cweight[ind])
    }
    rEstOld <- rEst
  }

  ## return values
  if(noBoot == TRUE)  {
    if(se == TRUE && length(unique(x)) > 1)
      ses <- les:::seFast(x, y, q)
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
  un <- !any(duplicated(pvalSort))

  if(un == TRUE)  {
    cdf <- cumsum(weightSort)
    cdf <- cdf/cdf[nPval]
    cdf <- cdf - (cdf-c(0, cdf[-nPval]))/2
  }
  else  {
    cdf <- les:::cdfDuplicates(pvalSort, weightSort)
  }
  res <- list(pval=pvalSort, cdf=cdf, unique=un)

  return(res)
}


##################################################
## cdfDuplicates
##################################################
cdfDuplicates <- function(pvalSort, weightSort)  {

  nPval <- length(pvalSort)
  tabCount <- rle(pvalSort)$lengths
  tabWeight <- table(pvalSort, weightSort)
  weightSum <- tabWeight %*% sort(unique(weightSort))
  cdf1 <- cumsum(weightSum)
  cdf1r <- rep.int(cdf1, tabCount)
  cdf2r <- rep.int(c(0, cdf1[-length(cdf1)]), tabCount)
  cdf2r <- cdf2r/cdf1r[nPval]
  cdf1r <- cdf1r/cdf1r[nPval]
  cdf <- cdf1r - (cdf1r - cdf2r)/2

  return(cdf)
}


##################################################
## grenanderPass
##################################################
grenanderPass <- function(x, y, un)  {

  if(un == TRUE)  {
    z <- GSRI:::grenanderInterp(x, y)
  }
  else  {
    ind <- !duplicated(x)
    z <- GSRI:::grenanderInterp(x[ind], y[ind])
    z <- rep.int(y[ind], rle(x)$lengths)
  }

  return(z)
}


##################################################
## seFast
##################################################
seFast <- function(x, y, b)  {
  
  si     <- sum((y-b*x)^2)/(length(x)-1)
  se     <- as.numeric(sqrt(si/(t(x)%*%x)))
  
  return(se)
}


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


##################################################
## gsri
##################################################
gsri <- function(pval, grenander=FALSE, se=TRUE, custom=FALSE)  {

  cweight <- rep(1, length(pval))
  res <- les:::fitGsri(pval, NULL, cweight, length(pval),
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


ind2log <- function(ind, n)  {

  log <- logical(n)
  log[ind] <- TRUE

  return(log)
}


##################################################
## slopeWeight
##################################################
slopeWeight <- function(x, y, c)  {
  s <- t(x) %*% c
  b <- ((s %*% y) / (s %*% x))[1]
  return(b)
}


qrSlope <- function(x, y, w)  {
  b <- as.numeric(stats:::lm.wfit(x, y, w)$coefficients)
  return(b)
}


diagSquare <- function(w, n)  {
  y <- rep(0, n*n)
  y[1+0:(n-1)*(n+1)] <- w
  dim(y) <- c(n, n)
  return(y)
}
