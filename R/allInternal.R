######################################################################
## internal functions for les package
######################################################################

######################################################################
## calcSingle
######################################################################
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
      res <- c(NA, NA, nValidProbes, NA)
    else
      res <- c(NA, NA)
  }
  
  return(res)
}


######################################################################
## fitGsri
######################################################################
fitGsri <- function(pval, index=NULL, cweight,
                    nValidProbes, grenander, se, custom)  { ## boot flag?

  noBoot <- is.null(index)
  if(noBoot == FALSE) {
    pval <- pval[index]
    cweight <- cweight[index]
  }
  cdf <- les:::wcdfGrenander(pval, cweight, nValidProbes, grenander, custom)
  res <- les:::itLinReg(cdf$pval, cdf$cdf, cweight, nValidProbes, se, custom, noBoot)
  if(noBoot == TRUE)
    res[4] <- cdf$l0

  return(res)
}


######################################################################
## itLinReg
######################################################################
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
    rEst <- max(c(rEst, 1))
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
    if(se == TRUE && length(unique(x)) > 1)  {
      xi <- x[ind]
      dim(xi) <- c(length(ind), 1)
      fit <- stats::lm.wfit(xi, y[ind], cweight[ind])
      ses <- les:::seSlopeWeight(fit, length(ind))
    }
    else  {
      ses <- NA
    }
    res <- c(1 - min(q, 1), ses, nValidProbes, NA)
  }
  else  {
    res <- 1 - min(q, 1, na.rm=TRUE)
  }
  return(res)
}


######################################################################
## wcdf2
######################################################################
wcdf2 <- function(pval, weight)  {
  
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


######################################################################
## wcdfGrenander
######################################################################
wcdfGrenander <- function(pval, cweight, nValidProbes, grenander, custom) {

  cdf <- les:::wcdf2(pval, cweight)
  cdf$l0 <- NA
  cdf$grenander <- grenander
  if(grenander == TRUE) {
    l0 <- les:::itLinReg(cdf$pval, cdf$cdf, cweight, nValidProbes, FALSE, custom, FALSE) ## noBoot=FALSE!
    cdf$cdf <- les:::modGrenander(cdf$pval, cdf$cdf, 1-l0, cdf$unique, FALSE)
    cdf$l0 <- l0
  }

  return(cdf)
}


######################################################################
## cdfDuplicates
######################################################################
cdfDuplicates <- function(pvalSort, weightSort)  {

  nPval <- length(pvalSort)
  #tabCount <- rle(pvalSort)$lengths ## this line seems buggy, 2014-03-25
  tabCount <- table(pvalSort)
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


######################################################################
## modGrenander
######################################################################
modGrenander <- function(x, y, eta0, unique=TRUE, check=TRUE) {
  
  ## for type="lcm" only, modification of "gcmlcm" (fdrtool)
  if(unique == FALSE) {
    ind <- !duplicated(x)
    x0 <- x
    x <- x[ind]
    y <- y[ind]
  }
  if(check == TRUE) {
    if(is.unsorted(x))
      stop("The x values must be arranged in sorted order.")
    if(any(duplicated(x)))
      stop("No duplicated x values allowed.")
  }
  yLower <- eta0*x
  yUpper <- 1 - eta0*(1 - x) ## == (1-eta0) + p*eta0
  indLower <- y < yLower
  indUpper <- y > yUpper
  y[indLower] <- yLower[indLower]
  y[indUpper] <- yUpper[indUpper]
  dx <- diff(x)
  dy <- diff(y)
  rawslope <- -dy/dx
  n <- length(rawslope)
  if(n == 1) {
    slope <- rawslope
  }
  else {
    slope <- .C("C_isomean", as.double(rawslope), as.double(dx),
                as.integer(n), ghat=double(n), PACKAGE="fdrtool", DUP=FALSE)$ghat
  }
  slope <- -slope
  yi <- y[1] + c(0, cumsum(dx*slope))

  if(unique == FALSE) {
    yi <- rep.int(yi, rle(x0)$lengths)
  }

  return(yi)
}


######################################################################
## seFast
######################################################################
seFast <- function(x, y, b)  {
  
  si     <- sum((y-b*x)^2)/(length(x)-1)
  se     <- as.numeric(sqrt(si/(t(x)%*%x)))
  
  return(se)
}


######################################################################
## mcapply
######################################################################
mcsapply <- function(X, FUN, ..., mc.cores=NULL, simplify=TRUE)  {

  mcLoaded <- "parallel" %in% .packages() &&
  any("mclapply" %in% objects("package:parallel"))

  if(mcLoaded == TRUE && !is.null(mc.cores))  {
    res <- parallel::mclapply(X, FUN, ...,  mc.cores=mc.cores)
    res <- sapply(res, c, simplify=simplify)
  }
  else  {
    res <- sapply(X, FUN, ..., simplify=simplify)
  }

  return(res)
}


######################################################################
## gsri
######################################################################
gsri <- function(pval, grenander=FALSE, se=TRUE, custom=FALSE)  {

  cweight <- rep(1, length(pval))
  res <- les:::fitGsri(pval, NULL, cweight, length(pval),
                       grenander=grenander, se=se, custom=custom)
  res[4] <- res[1]*res[3]
  names(res) <- c("GSRI", "se", "n", "nReg")

  return(res)
}


######################################################################
## log2ind
######################################################################
log2ind <- function(log)  {

  ind <- seq(along=log)[log]
  ind <- as.integer(ind)

  return(ind)
}


######################################################################
## ind2log
######################################################################
ind2log <- function(ind, n)  {

  log <- logical(n)
  log[ind] <- TRUE

  return(log)
}


######################################################################
## reg2log
######################################################################
reg2log <- function(reg, pos, chr)  {

  if(length(reg$chr) == 0)
    indChr <- !logical(length(pos))
  else
    indChr <- chr %in% reg$chr
  ind <- indChr & pos >= reg$start & pos <= reg$end

  return(ind)
}


######################################################################
## slopeWeight
######################################################################
slopeWeight <- function(x, y, c)  {
  
  s <- t(x) %*% c
  b <- ((s %*% y) / (s %*% x))[1]
  
  return(b)
}


######################################################################
## qrSlope
######################################################################
qrSlope <- function(x, y, w)  {
  
  b <- as.numeric(stats:::lm.wfit(x, y, w)$coefficients)
  
  return(b)
}


######################################################################
## seSlopeWeight
######################################################################
seSlopeWeight <- function(fit, n=length(fit$weights))  {
  
  R <- fit$qr$qr[1]^-2  ## = chol2inv(fit$qr$qr[1])
  rss <- sum(fit$weights*fit$residuals^2)
  resvar <- rss/(n - 1)
  se <- sqrt(R*resvar)
  
  return(se)
}


######################################################################
## diagSquare
######################################################################
diagSquare <- function(w, n)  {
  
  y <- rep(0, n*n)
  y[1+0:(n-1)*(n+1)] <- w
  dim(y) <- c(n, n)
  
  return(y)
}

######################################################################
## scaleNorm
######################################################################
scaleNorm <- function(x)  {

  x <- x - min(x, na.rm=TRUE)
  m <- max(x, na.rm=TRUE)
  if(m != 0)
    x <- x/m
  
  return(x)
}


######################################################################
## xvalWeight
######################################################################
xvalWeight <- function(distance, win)  {

  weighting <- getOption("weighting")
  weight <- weighting(distance, win)
  weight[distance == 0] <- 0
  
  return(weight)
}


######################################################################
## optimalSingleRegion
######################################################################
optimalSingleRegion <- function(i, reg, object, winSize, fdr, method,
                                scaling, nCores, verbose, ...)  {
  
  nWin <- length(winSize)
  chi2 <- vector("numeric", nWin)
  ind <- les:::reg2log(reg[i, ], object@pos, object@chr)
  resc <- les::Les(object@pos[ind], object@pval[ind])

  if(verbose == TRUE)
    print(sprintf("%s %d/%d", "Region", i, nrow(reg)))
  if(sum(ind) == 0)
    warning(sprintf("%s %d", "No probes in region", i))
  
  for(w in 1:nWin)  {
    lesi <- les::estimate(resc, winSize[w], weighting=les:::xvalWeight,
                          grenander=object@grenander, se=FALSE,
                          minProbes=object@minProbes, method=method,
                          nCores=nCores, verbose=FALSE)@lambda
    
    ## take only non-NAs into account
    indValid <- !is.na(lesi) & !is.na(fdr[ind])
    
    chi2[w] <- sum((scaling(lesi[indValid]) -
                    scaling(1 - fdr[ind][indValid]))^2)/sum(indValid)
  }
  
  return(chi2)
}


######################################################################
## setState
######################################################################
setState <- function(state, flag)  {

  if(!(flag %in% state))
    state <- c(state, flag)

  return(state)
}


######################################################################
## checkState
######################################################################
checkState <- function(state, flag, text=as.character(flag))  {

  status <- flag %in% state
  if(all(!status))
    stop(sprintf("'%s' %s", text, "must be called first."))

  return(status)
}


getArgs <- function(name, first=NULL, last=NULL, ...) {

  ind <- which(names(...) %in% name)[1] ## [[]] allows only one element for indexing
  middle <- if(length(ind) != 0) ...[[ind]] else NULL
  args <- c(first, middle, last)
  args <- args[!duplicated(names(args))]
  
  return(args)
}

