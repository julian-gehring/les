##################################################
## optimalKernel
##################################################
#setMethod("optimalKernel", "Les",
#          function(object, winSize, regions, verbose=FALSE)  {

## name already used by default package
optimize <- function(object, winSize, regions, fdr, method, verbose=FALSE, nCores=NULL)  {

  if(missing(regions))  {
    if(length(object@regions) != 0)
      regions <- object@regions
    else
      stop("'regions' must be specified.")
  }

  nReg <- nrow(regions)
  nWin <- length(winSize)

  ## get fdr
  if(is.character(fdr))  {
    fdrMethod <- c("lfdr", "qval")
    choice <- pmatch(fdr, fdrMethod)
    if(!is.na(choice))
      fdr <- fdrtool::fdrtool(object@pval, "pvalue", verbose=FALSE, plot=FALSE)[[fdrMethod[choice]]]
    else
      stop("'fdr' is not a valid fdr method.")
  }
  else  {
    if(length(fdr) != length(object@pos))
      stop("'fdr' must have the same size as 'pos'.")
  }

  ## set weighting in options
  if(length(object@weighting) == 0)
    stop("'weighting' function must be set in the 'estimate' method.")
  oldOpt <- getOption("weighting")
  options(weighting=object@weighting)

  ## create data object
  if(missing(method))
    method <- object@method

  chi2 <- sapply(1:nReg, optimalSingleRegion, regions, object, winSize, fdr, method, nCores, verbose=verbose)

  rownames(chi2) <- winSize
  colnames(chi2) <- rownames(regions)

  ## unset weighting in options
  options(weighting=NULL)

  return(chi2)
}



optimalKernel <- function(object, winSize, regions, fdr, verbose=FALSE, scaling=scaleNorm)  {
  
  if(missing(regions))
    regions <- rep(TRUE, length(object@pos))
  nReg <- 1
  if(missing(fdr))
    fdr <- fdrtool::fdrtool(object@pval, "pvalue", verbose=FALSE, plot=FALSE)$lfdr

  if(length(object@weighting) == 0)
    stop("'weighting' function must be set by the 'estimate' method.")
  options(weighting=object@weighting)

  chi2 <- matrix(NA, nReg, length(winSize))
  ##for(r in 1:nReg)  {
  ##ind <- inVector(object@pos, regions[r,1], regions[r,2])
  ind <- regions
  #suppressWarnings(fdr <- fdrtool::fdrtool(object@pval[ind], "pvalue",
                          #verbose=FALSE, plot=FALSE))
  ## if(length(unique(object@chr[ind])) > 1)
  ## stop("regions covers more than 1 chromosome.")
  resc <- create(object@pos[ind], object@pval[ind])
  for(w in 1:length(winSize))  {
    if(verbose == TRUE)
      print(round(w/length(winSize), 2))
    lesi <- estimate(resc, winSize[w], weighting=xvalWeight, nCores=2, method="qr")["lambda"]
    indValid <- !is.na(lesi) & !is.na(fdr[ind])
    chi2[1,w] <- sum((scaling(lesi[indValid]) - scaling(1 - fdr[ind][indValid]))^2)/sum(indValid)
    ##browser()
    ##}
  }
  colnames(chi2) <- winSize

  ## unset weighting in options
  options(weighting=NULL)
  
  return(chi2)
}


xvalWeight <- function(distance, win)  {

  weighting <- getOption("weighting")
  weight <- weighting(distance, win) ## problematic ## NOT GOOD !!! ##
  weight[distance == 0] <- 0
  
  return(weight)
}


scaleNorm <- function(x)  {

  x <- x - min(x, na.rm=TRUE)
  if(max(x, na.rm=TRUE) == 0)
    x <- x + 1e-14
  x <- x/max(x, na.rm=TRUE)

  return(x)
}


reg2log <- function(reg, pos, chr)  {

  if(length(reg$chr) != 0)
    indChr <- chr %in% reg$chr
  else
    indChr <- !logical(length(pos))
  ind <- pos[indChr] >= reg$start & pos[indChr] <= reg$end

  return(ind)
}


optimalSingleRegion <- function(i, reg, object, winSize, fdr, method, nCores, verbose, ...)  {
  
  nWin <- length(winSize)
  chi2 <- vector("numeric", nWin)
  ind <- reg2log(reg[i, ], object@pos, object@chr)
  resc <- create(object@pos[ind], object@pval[ind])  ## chr not needed since only on 1?

  if(verbose == TRUE)
    print(sprintf("%s %d/%d", "Region", i, nrow(reg)))
  
  for(w in 1:nWin)  {
    lesi <- estimate(resc, winSize[w], weighting=les:::xvalWeight,
                     grenander=object@grenander, se=FALSE,
                     minProbes=object@minProbes, method=method,
                     nCores=nCores, verbose=FALSE)@lambda
    
    ## take only non-NAs into account
    indValid <- !is.na(lesi) & !is.na(fdr[ind])  ## ind !
    
    chi2[w] <- sum((les:::scaling(lesi[indValid]) -
                    les:::scaling(1 - fdr[ind][indValid]))^2)/sum(indValid)
  }
  
  return(chi2)
}
