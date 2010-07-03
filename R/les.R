##################################################
## create
##################################################
create <- function(pos, pval, chr)  {

  ## default values
  if(missing(chr))
    chr <- rep(0, length(pos))

  ## check inputs
  if(length(pos) != length(pval))
    stop("'pos' and 'pval' must have the same length.")
  if(length(pos) != length(chr))
    stop("'pos' and 'chr' must have the same length.")
  if(any(is.na(pos)))
    stop("'pos' must not contain NAs.")
  if(any(pos %% 1 != 0))
    stop("'pos' must be a vector of integers.")

  ## throw out NAs in pval
  indValid <- !is.na(pval)
  pos <- as.integer(pos[indValid])
  pval <- round(pval[indValid], 12)
  chr <- factor(chr[indValid])

  ## sort
  ord <- order(chr, pos)
  pos <- pos[ord]
  pval <- pval[ord]
  chr <- chr[ord]

  ## check for duplicated pos for each chr
  if(any(grepl("TRUE", rownames(table(duplicated(pos), chr)))))
    warning("'pos' contains duplicates on the same chr.")

  object <- new(Class="Les",
                pos=pos, pval=pval, chr=chr, nChr=nlevels(chr))
  
  return(object)
}
## ok ##


##################################################
## estimate
##################################################
setMethod("estimate", "Les",
          function(object, win, weighting=triangWeight,
                   grenander=FALSE, se=TRUE, nCores=NULL)  {
            
  ## check input
  if(class(object) != "Les")
    stop("'object' must be of class 'Les'")
  if(win %% 1 != 0)
    warning("'win' was rounded to nearest integer.")
  win <- as.integer(win)
  chrLevel <- levels(object@chr)

  ## for each chr
  for(c in 1:object@nChr)  {
    indChr <- object@chr == chrLevel[c]
    indProbes <- seq(win+1, win+length(object@pos[indChr]))
    pos <- c(rep(-Inf, win), object@pos[indChr], rep(Inf, win))
    pval <- c(rep(NA, win), object@pval[indChr], rep(NA, win))

    ## apply
    cs <- mcsapply(indProbes, calcSingle, pos, pval,
                   win, weighting, grenander, se, nBoot=FALSE, mc.cores=nCores)
    
    ## extract result
    object@lambda[indChr] <- cs[1, ]
    object@nProbes[indChr] <- as.integer(cs[3, ])
    if(se == TRUE)
      object@se[indChr] <- cs[2, ]
  }

  object@win <- win
  object@weighting <- weighting
  object@grenander <- grenander
  
  return(object)
}
)

## ok ##


##################################################
## calcSingle
##################################################
calcSingle <- function(ind0, pos, pval, win,
                       weighting, grenander, se,
                       nBoot, conf)  {

  ## cut down
  pos0 <- pos[ind0]
  indCut <- seq(ind0-win, ind0+win)
  posCut <- pos[indCut]
  pvalCut <- pval[indCut]
  
  distance <- abs(posCut - pos0)
  indValid <- distance <= win
  nValidProbes <- length(unique(pvalCut[indValid]))

  if(nValidProbes>1)  {
    ## if enough probes
    dis <- distance[indValid]
    weight <- weighting(dis, win)

    ## apply
    if(nBoot == FALSE)  {
      res <- fitGSRI(pvalCut[indValid], index=NULL, weight,
                     nValidProbes, grenander, se)
    }
    else  {
      bo <- boot::boot(pvalCut[indValid], fitGSRI, nBoot,
                 cweight=weight, nValidProbes=nValidProbes,
                 grenander=grenander, se=se)
      ## try not needed?
      ci <- try(boot::boot.ci(bo, conf, type="perc"), silent=FALSE)
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
## ok ## nValidProbes as input?


##################################################
## fitGSRI
##################################################
fitGSRI <- function(pval, index=NULL, cweight, nValidProbes, grenander, se)  {
  
  maxIter <- nValidProbes
  restOld <- 0
  rest <- restOld
  q <- 1
  noBoot <- is.null(index)
  if(noBoot == TRUE)
    cdf <- wcdf(pval, cweight, grenander)
  else
    cdf <- wcdf(pval[index], cweight[index], grenander)

  ## iterative fitting
  for(i in 1:maxIter)  {
    rest <- nValidProbes - ceiling(q*nValidProbes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(c(nValidProbes-1, rest))
    if(is.na(rest) || restOld == rest)
      break
    x <- cdf$pval[rest:nValidProbes] - 1
    y <- cdf$cdf[rest:nValidProbes] - 1
    if(length(unique(x)) == 1) # only for boot?
      break
    q <- GSRI:::slopeFast(x, y)
    restOld <- rest
  }

  ## return values
  if(noBoot == TRUE)  {
    if(se == TRUE)
      se <- seFast(x, y, q)
    else
      se <- NA
    res <- c(1 - min(q, 1), se, nValidProbes)
  }
  else  {
    res <- 1 - min(q, 1, na.rm=TRUE)
  }
  
  return(res)
}
## ok ## needs testing, nValidProbes


##################################################
## wcdf
##################################################
wcdf <- function(pval, weight, grenander)  {
  
  ord <- sort.list(pval, method="quick", na.last=NA)
  pvalSort <- pval[ord]
  weightSort <- weight[ord]

  indUnique <- !duplicated(pvalSort)
  uniquePval <- pvalSort[indUnique]
  uniqueWeight <- weightSort[indUnique]

  nUnique <- length(uniquePval)
  nProbes <- length(pvalSort)

  cdf <- cumsum(uniqueWeight)
  cdf <- cdf/cdf[nUnique]
    
  if(nProbes != nUnique)
    cdf <- rep.int(cdf, table(pvalSort))
  if(grenander == TRUE)  {
    if(nProbes != nUnique)  {
      jit <- round(runif(nProbes-nUnique, 2e-12, 5e-12), 12)
      pvalSort[!indUnique] <- pvalSort[!indUnique] + jit
    }
    cdf <- GSRI:::grenanderInterp(pvalSort, cdf)
  }
  cdf <- cdf - 0.5/nProbes  ## where to put this ??
  res <- list(pval=pvalSort, cdf=cdf)

  return(res)
}


##################################################
## seFast
##################################################
seFast <- function(x, y, b) {
  
  si     <- sum((y-b*x)^2)/(length(x)-1)
  se     <- as.numeric(sqrt(si/(t(x)%*%x)))
  
  return(se)
}
## ok ##


##################################################
## triangWeight
##################################################
triangWeight <- function(distance, win)  {
  
  weight <- 1 - distance/win  ## not normed
  ## does not matter, normed in cdf
    
  return(weight)
}
## ok ##


##################################################
## rectangWeight
##################################################
gaussWeight <- function(distance, win)  {

  weight <- dnorm(distance, sd=win/2)

  return(weight)
}
## ok ##


##################################################
## gaussWeight
##################################################
rectangWeight <- function(distance, win)  {

  n <- length(distance)
  weight <- rep(1/n, n)

  return(weight)
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
## ci
##################################################
setMethod("ci", "Les",
          function(object, subset, nBoot=100, conf=0.95, nCores=NULL)  {

  if(class(object) != "Les")
    stop("'object' must be of class 'Les'")
  if(missing(subset))
    subset <- rep(TRUE, length(object@pos))
    #subset <- seq(along=object@pos)
  if(!is.logical(subset) & is.vector(subset))
    subset <- ind2log(subset, length(object@pos))
  #if(is.matrix(subset))
  #  subset <- reg2log(subset, object@pos)
  nBoot <- as.integer(nBoot)
  win <- object@win
  chrLevel <- levels(object@chr)
  nChr <- length(chrLevel)

  bc <- c()  ## make better
  for(c in 1:nChr)  {
    indChr <- object@chr == chrLevel[c]
    indProbes <- seq(win+1, win+length(object@pos[indChr]))[subset[indChr]]
    pos <- c(rep(-Inf, win), object@pos[indChr], rep(Inf, win))
    pval <- c(rep(NA, win), object@pval[indChr], rep(NA, win))
    
    bs <- mcsapply(indProbes, calcSingle,pos, pval,
                   win, object@weighting, object@grenander, se=FALSE, nBoot, conf,
                   mc.cores=nCores)
    bc <- cbind(bc, bs)
  }
  
  object@ci <- data.frame(lower=bc[1, ], upper=bc[2, ])
  object@subset <- log2ind(subset)
  object@nBoot <- nBoot
  object@conf <- conf
  
  return(object)
}
)


##################################################
## regions
##################################################
setMethod("regions", "Les",
          function(object, limit=NULL, minLength=10, maxGap=100, verbose=FALSE)  {

  if(class(object) != "Les")
    stop("'object' must be of class 'Les'")
  if(is.null(limit))  {
    if(length(object@theta) == 0)
      stop("'limit' must be specified")
    else
      limit <- object@theta
  }
  if(is.na(limit) || limit < 0 || limit > 1)
    stop("'limit' must be in the range [0,1].")
  lambda <- object@lambda
  indNa <- is.na(lambda)
  lambda[indNa] <- 0
  ind <- c(0, lambda >= limit, 0)
  indMax <- which(diff(object@pos) > maxGap)
  ind[indMax] <- 0
  ind[indMax+1] <- 0
  d <- diff(ind)
  begin <- which(d == 1)
  end <- which(d == -1) - 1
  nProbes <- end - begin +1
  
  indSame <- (object@chr[begin] == object@chr[end]) & (nProbes >= minLength)
  begin <- begin[indSame]
  end <- end[indSame]
  nProbes <- nProbes[indSame]

  if(length(begin) != length(end))
    stop("Unequal numbers of region bounderies")
  if(verbose == TRUE)
    print(sprintf("%d %s%g", length(begin), "regions found with lambda>=", limit))

  ri <- matrix(NA, length(begin), 2)
  for(i in 1:length(begin))  {
    ri[i, ] <- gsri(object@pval[begin[i]:end[i]], se=TRUE)[c(1,2)]
  }
  rs <- ri[ ,1]/ri[ ,2]
  
  size <- as.integer(object@pos[end]-object@pos[begin]+1)
  regions <- data.frame(chr=factor(object@chr[begin]),
                        start=object@pos[begin], end=object@pos[end],
                        size=size, nProbes=nProbes, ri=round(ri[ ,1], 4),
                        se=round(ri[ ,2], 4), rs=round(rs, 4))
  regions <- regions[order(rs, ri[ ,1], nProbes, size, decreasing=TRUE), ]

  object@regions <- regions
  object@limit <- limit
  object@minLength <- as.integer(minLength)
  object@maxGap <- as.integer(maxGap)
  
  return(object)
}
)


##################################################
## gsri
##################################################
gsri <- function(pval, grenander=FALSE, se=TRUE)  {

  cweight <- rep(1, length(pval))
  res <- fitGSRI(pval, NULL, cweight, length(pval), grenander, se)
  res <- c(res, res[1]*res[3])
  names(res) <- c("GSRI", "se", "n", "nReg")

  return(res)
}


##################################################
## cutoff
##################################################
setMethod("cutoff", "Les",
          function(object, grenander=FALSE, verbose=FALSE)  {
  
  if(class(object) != "Les")
    stop("'object' must be of class 'Les'")
  
  pval <- object@pval
  nProbes <- length(pval)
  erg <- gsri(pval, grenander, se=FALSE)
  nSigProbes <- erg[1]*nProbes
  nSigLower <- ceiling(nSigProbes)
  cutoff <- c(NA, sort(object@lambda, decreasing=TRUE))[nSigLower+1]
  ## +1: access right probe due to ceiling !
  
  if(verbose == TRUE)
    print(sprintf("%g %s%g", nSigLower,
                  "significant probes estimated with limit lambda>=", cutoff))
  
  object@nSigProbes <- nSigProbes ## round?
  object@theta <- cutoff
  
  return(object)
}
)


##################################################
## plot
##################################################
setMethod("plot", "Les",
          function(x, y, chr, region=FALSE, limit=TRUE,
                   xlim, ylim=c(0,1),
                   error="none", semSpread=1.96,
                   patchCol="lightgray", borderCol="black",
                   markerCol="red", regionCol,
                   regionLty=2, regionLw=3, regionCex=NULL,
                   probeCol="black", probeCex=0.5, probePch=20)  {

  pos <- x@pos
  lambda <- x@lambda
  
  if(error == "ci" && length(x@ci) == 0)  {
    warning("CI not computed so far")
    error <- "none"
  }

  if(missing(chr))  {
    if(x@nChr != 1)
      stop("Please specify 'chr'.")
    else
      indChr <- rep(TRUE, length(x@pos))
  }
  else  {
    indChr <- x@chr %in% chr
  }
  
  if(!missing(y))  {
    par(mar=c(2.1,4.1,4.1,2.1))
    par(fig=c(0,1,0,0.3))
    plot.new()
    plot.window(x=xlim, y=c(0, 2))
    abline(h=0.5, col="gray32")
    indAnno <- reg2ind(y$reg, xlim)
    for(i in indAnno)  {
      polygon(c(y$pos[ ,i], rev(y$pos[ ,i])),
              c(0, 0, 1, 1), col=y$col[i])
    }
    text(x=0.5*(anno$pos[2,indAnno]+anno$pos[1,indAnno]),
         y=0.5,
         label=anno$name[indAnno],
         adj=c(0.5,0.5)
         )
    text(x=xlim[2], y=0.55,
         label="DNA", adj=c(0.5,0))
    ## for next plot
    par(fig=c(0,1,0.2,1), new=TRUE)
  }
  
  if(missing(xlim))
    xlim <- range(pos)
  ind <- (pos >= xlim[1]) & (pos <= xlim[2])  ## needed for plotting?
  
  if(is.logical(limit) && limit == TRUE)  {
    if(length(x@theta) != 0)
      limit <- x@theta
    else  {
                                        #warning("'cutoff' not estimated")
      limit <- FALSE  ## skip warning??
    }
  }
  if(region == TRUE && length(x@regions) != 0)
    regions <- x@regions
  if(!is.logical(region))
    regions <- region
  
  sig <- lambda >= limit
  plot(pos, lambda, type="n", xlim, ylim,
       xlab="Probe position", ylab=expression(Lambda))
  if(is.numeric(limit))
    abline(h=limit, col="gray")
  abline(h=0, col="lightgray")
  sfrac <- min(min(diff(pos[ind]))/(xlim[2]-xlim[1])/2, 0.01)
  suppressWarnings(
    switch(match(error, c("ci", "se", "none")),  {
      ss <- x@subset
      ci <- x@ci
      plotCI(pos[ss], lambda[ss],
             ui=ci$upper, li=ci$lower,
             gap=0, pch=".", col="azure4",
             add=TRUE, sfrac=sfrac)
    },  {
      se <- x@se
      plotCI(pos, lambda,
             uiw=se*semSpread, liw=se*semSpread,
             gap=0, pch=".", col="azure4",
             add=TRUE, sfrac=sfrac)
    }
           )
                   )
  
  points(pos, lambda, type="o",
         pch=probePch, col=probeCol, cex=probeCex)
  if(is.numeric(limit))  {
    points(pos[sig], lambda[sig],
           pch=probePch, col=markerCol, cex=probeCex)
  }
  
  if(!is.logical(region) && !is.numeric(region))
    stop("'region' must be logical or numeric")
  if(is.logical(region) && region == TRUE)
    region <- x@regions
  if(!is.logical(region) || region == TRUE)  {
    regionCol <- rep(brewer.pal(8, "Set1"), each=2)
    abline(v=region, col=regionCol)
  }
}
          )


##################################################
## [
##################################################
setMethod("[", "Les",
          function(x, i, j, drop)  {
            value <- slot(x, i)
            return(value)
          }
)
## ok ##


##################################################
## <-
##################################################
setReplaceMethod("[", "Les",
                 function(x, i, j, value)  {
                   slot(x, i) <- value
                   validObject(x)
                   return(x)
                 }
)
## ok ##  is this function needed by the user ? ##


##################################################
## show
##################################################
setMethod("show", "Les",
          function(object)  {
            
  pos <- object@pos
  lambda <- object@lambda
  nChr <- object@nChr
  cat("** Object of class 'Les' **\n")
  cat(sprintf("* %d %s %d %s\n", length(pos),
              "probes on", nChr, "chromosomes"))
  if(length(lambda) != 0)
    cat(sprintf("* %s [%g, %g] %s %d\n", "Lambda in range",
                min(lambda, na.rm=TRUE), max(lambda, na.rm=TRUE),
                "with kernel size", object@win))
  if(length(object@ci) != 0)
    cat(sprintf("* %g %s %d %s\n", object@conf, "confidence intervals computed for",
                length(object@subset), "probes"))
  if(length(object@nSigProbes) != 0)
    cat(sprintf("* %d %s %g\n", ceiling(object@nSigProbes),
                "regulated probes estimated for lambda >=",
                object@theta))
  if(length(object@regions) != 0)
    cat(sprintf("* %d %s\n", nrow(object@regions), "regions detected"))
}
          )


##################################################
## summary
##################################################
setMethod("summary", "Les",
          function(object, ...)  {
            show(object, ...)
            cat("* Available slots:\n")
            sn <- slotNames(object)
            ind <- sapply(sn, function(sn) length(object[sn])) > 0
            print(sn[ind])
          }
          )


##################################################
## optimalKernel
##################################################
setMethod("optimalKernel", "Les",
          function(object, winSize, regions, verbose=FALSE)  {
  
  if(missing(regions))
    regions <- rep(TRUE, length(object@pval))
  #if(is.matrix(regions))
  #  nReg <- nrow(regions)
  #else
    nReg <- 1
  #if(ncol(regions) != 2)
  #  stop("'regions' must have 2 columns")

  if(length(object@weighting) == 0)
    stop("'weighting' function must be set by the 'estimate' method")
  options(weighting=object@weighting)
  
  ## check if existing

  chi2 <- matrix(NA, nReg, length(winSize))
  for(r in 1:nReg)  {
    #ind <- inVector(object@pos, regions[r,1], regions[r,2])
    ind <- regions
    fdr <- fdrtool::fdrtool(object@pval[ind], "pvalue",
                            verbose=FALSE, plot=FALSE, cutoff="pct0")
    if(length(unique(object@chr[ind])) > 1)
      stop("regions covers more than 1 chromosome.")
    resc <- create(object@pos[ind], object@pval[ind])
    for(w in 1:length(winSize))  {
      lesi <- estimate(resc, winSize[w], weighting=xvalWeight, nCores=FALSE)["lambda"]
      #chi2[r,w] <- sum((scaling(lesi) + scaling(fdr$qval) - 1)^2)
      chi2[r,w] <- sum((scaling(lesi) - scaling(1 - fdr$lfdr))^2)
    }
  }
  colnames(chi2) <- winSize
  
  return(chi2)
}
)


##################################################
## xvalWeight
##################################################
xvalWeight <- function(distance, win)  {

  weighting <- getOption("weighting")
  weight <- weighting(distance, win) ## problematic ## NOT GOOD !!! ##
  weight[distance == 0] <- 0
  
  return(weight)
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


reg2log <- function(reg, x)  {

  if(ncol(reg) != 2)
    stop("'reg' must have 2 columns")
  log <- logical(length(x))
  for(i in 1:nrow(reg))  {
    ind <- x >= min(reg[i, ]) & x <= max(reg[i, ])
    log[ind] <- TRUE
  }

  return(log)
}
## ok ##


reg2ind <- function(reg, x)  {

  log <- reg2log(reg, x)
  ind <- log2ind(log)

  return(ind)
}
## ok ##


log2reg <- function(log)  {

  ind <- c(FALSE, log, FALSE)
  d <- diff(ind)
  s <- seq(along=d)
  begin <- s[d == 1]
  end <- s[d == -1] - 1
  #ind <- begin != end
  reg <- cbind(begin, end)
  colnames(reg) <- c("begin", "end")

  return(reg)
}
## ok ##


ind2reg <- function(ind, n)  {

  log <- ind2log(ind, n)
  reg <- log2reg(log)

  return(reg)
}
## ok ##


scaling <- function(x)  {

  x <- x - min(x)
  x <- x/max(x)

  return(x)
}
## ok ##


inVector <- function(pos, start, end)  {
  if(length(start) != length(end))
    stop("not same length")
  sig <- logical(length(pos))
  for(i in 1:length(start))  {
    ind <- pos >= start[i] & pos <= end[i]
    sig[ind] <- TRUE
  }
  return(sig)
}
## ok ##

setMethod("exportLambda", "Les",
          function(object, chr, file, range,
                   description="Lambda", precision=4)  {

  if(class(object) != "Les")
    stop("'object' must be of class 'Les'")
            
  header <- c(sprintf("%s%s%s","track name=\"", description,
                      "\" type=wiggle_0 viewLimits=0:1 autoScale=off"),
              sprintf("%s%s %s", "variableStep chrom=", chr,  "span=1"))

  if(missing(range))  {
    ind <- object@chr == chr & !is.na(object@lambda)
  }
  else  {
    ind <- object@chr == chr & object@pos <= min(range) &
    object@pos >= max(range) & !is.na(object@lambda)
  }

  values <- format(round(object@lambda[ind], precision), scientific=16)
  df <- data.frame(pos=object@pos[ind], value=values)
  # not good, why this?
  indDup <- duplicated(df$pos)
  df <- df[!indDup, ]
  indValid <- !is.na(df$value)
  df <- df[indValid, ]

  write(header, file)
  write.table(df, file, append=TRUE, quote=FALSE, sep="\t",
              row.names=FALSE, col.names=FALSE)

}
)
