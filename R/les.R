##################################################
## create
##################################################
create <- function(pos, pval, chr)  {

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
  #if(!is.numeric(pval) || min(pval, na.rm=TRUE) < 0 || max(pval, na.rm=TRUE) > 1)
  #  stop("'pval' must contain numerics in the range [0,1].")
  ## bring back later !!!!

  ## throw out NAs in pval
  indValid <- !is.na(pval)
  pos <- as.integer(pos[indValid])
  # pval <- round(pval[indValid], 12)
  chr <- factor(chr[indValid])

  ## sort
  ord <- order(chr, pos)
  pos <- pos[ord]
  pval <- pval[ord]
  chr <- chr[ord]

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
                   se=FALSE, minProbes=3, method="la", nCores=NULL, ...)  {
            
  ## check input
  if(missing(win))
    stop("'win' must be specified.")
  if(win < 1)
    stop("'win' must be a positive integer.")
  if(win %% 1 != 0)
    warning("'win' was rounded to nearest integer.")
  custom <- pmatch(method, c("la", "qr"), NA) == 1
  if(is.na(custom))
    stop("'method' must be 'la' or 'qr'")
  win <- as.integer(win)
  minProbes <- as.integer(minProbes - 1)
  chrLevel <- levels(object@chr)

  ## silent default for grenander
  arg <- list(...)
  grenander <- FALSE
  if(any(names(arg) == "grenander"))
    grenander <- arg[[which(names(arg) == "grenander")]]

  ## for each chr
  for(c in 1:object@nChr)  {
    indChr <- object@chr == chrLevel[c]
    indProbes <- seq(win+1, win+length(object@pos[indChr]))
    pos <- c(rep(-Inf, win), object@pos[indChr], rep(Inf, win))
    pval <- c(rep(NA, win), object@pval[indChr], rep(NA, win))

    ## apply
    cs <- mcsapply(indProbes, calcSingle, pos, pval,
                   win, weighting, grenander, se, nBoot=FALSE,
                   minProbes=minProbes, custom=custom, mc.cores=nCores)
    
    ## extract result
    object@lambda[indChr] <- cs[1, ]
    object@nProbes[indChr] <- as.integer(cs[3, ])
    if(se == TRUE)
      object@se[indChr] <- cs[2, ]
  }

  object@win <- win
  object@weighting <- weighting
  object@grenander <- grenander
  object@minProbes <- minProbes
  object@method <- method
  
  return(object)
}
)

## ok ##


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
      res <- fitGSRI(pvalCut[indValid][indWeight], index=NULL,
                     weight[indWeight], nWeight, grenander, se,
                     custom)
    }
    else  {
      bo <- boot::boot(pvalCut[indValid][indWeight], fitGSRI, nBoot,
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
## fitGSRI
##################################################
fitGSRI <- function(pval, index=NULL, cweight,
                    nValidProbes, grenander, se, custom)  {
  
  maxIter <- nValidProbes
  restOld <- 0
  rest <- restOld
  q <- 1
  noBoot <- is.null(index)
  if(noBoot == TRUE)
    cdf <- wcdf2(pval, cweight, grenander)
  else
    cdf <- wcdf2(pval[index], cweight[index], grenander)
  if(custom == TRUE)
    cw <- diagSquare(cweight, nValidProbes)
  if(any(cdf$cdf < 0))
    stop("weights < 0")
  x <- cdf$pval - 1
  y <- cdf$cdf - 1
  ## iterative fitting
  for(i in 1:maxIter)  {
    rest <- nValidProbes - ceiling(q*nValidProbes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(c(nValidProbes-1, rest))
    if(is.na(rest) || restOld == rest)
      break
#    if(length(unique(x)) == 1) # only for boot?
#      break
    ind <- rest:nValidProbes
    nRest <- length(ind)
    if(custom == TRUE)  {
      q <- slopeWeight(x[ind], y[ind], cw[ind,ind])
    }
    else  {
      xi <- x[ind]
      dim(xi) <- c(nRest, 1)
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
## ok ## needs testing, nValidProbes


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
## triangWeight
##################################################
triangWeight <- function(distance, win)  {
  
  weight <- 1 - abs(distance)/win
    
  return(weight)
}
## ok ##


##################################################
## gaussWeight
##################################################
gaussWeight <- function(distance, win)  {

  weight <- stats::dnorm(distance, sd=win/2)

  return(weight)
}
## ok ##


##################################################
## rectangWeight
##################################################
rectangWeight <- function(distance, win)  {

  n <- length(distance)
  weight <- rep(1/n, n)

  return(weight)
}
## ok ##


##################################################
## epWeight
##################################################
epWeight <- function(distance, win)  {
  
  weight <- 0.75*(1 - (distance/win)^2)
    
  return(weight)
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
## ok ##


##################################################
## ci
##################################################
setMethod("ci", "Les",
          function(object, subset, nBoot=100, conf=0.95,
                   nCores=NULL, ...)  {

  if(missing(subset))
    subset <- rep(TRUE, length(object@pos))
    #subset <- seq(along=object@pos)
  if(!is.logical(subset) & is.vector(subset))
    subset <- ind2log(subset, length(object@pos))
  #if(is.matrix(subset))  ## TO DO: implement
  #  subset <- reg2log(subset, object@pos)
  nBoot <- as.integer(nBoot)
  win <- object@win
  custom <- object@method == 1
  chrLevel <- levels(factor(object@chr[subset]))
  nChr <- length(chrLevel)

  indRes <- seq(along=subset)[subset]
  bc <- c()
  for(c in 1:nChr)  {
    indChr <- object@chr == chrLevel[c]
    indProbes <- seq(win+1, win+length(object@pos[indChr]))[subset[indChr]]
    pos <- c(rep(-Inf, win), object@pos[indChr], rep(Inf, win))
    pval <- c(rep(NA, win), object@pval[indChr], rep(NA, win))
    bs <- mcsapply(indProbes, calcSingle, pos, pval,
                   win, object@weighting, object@grenander, se=FALSE,
                   nBoot, conf, minProbes=object@minProbes, custom=custom,
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
          function(object, limit=NULL, minLength=2,
                   maxGap=Inf, verbose=FALSE, ...)  {

  if(is.null(limit))  {
    if(length(object@theta) == 0)
      stop("'limit' must be specified.")
    else
      limit <- object@theta
  }
  else  {
    if(is.na(limit) || limit <= 0 || limit >= 1)
      stop("'limit' must be in the range [0,1].")
  }
  lambda <- object@lambda
  indNa <- is.na(lambda)
  lambda[indNa] <- 0
  ind <- c(0, lambda >= limit, 0)
  indMax <- which(diff(object@pos) > maxGap)
  ind[c(indMax, indMax+1)] <- 0
  d <- diff(ind)
  begin <- which(d == 1)
  end <- which(d == -1) - 1
  nProbes <- end - begin + 1
  
  indSame <- (object@chr[begin] == object@chr[end]) & (nProbes >= minLength)
  begin <- begin[indSame]
  end <- end[indSame]
  nProbes <- nProbes[indSame]

  if(length(begin) != length(end))
    stop("Unequal numbers of region bounderies found")
  if(verbose == TRUE)
    print(sprintf("%d %s%g", length(begin), "regions found for lambda>=", limit))

  ri <- matrix(NA, length(begin), 2)
  if(length(begin) > 0)  {
    for(i in 1:length(begin))  {
      ri[i, ] <- gsri(object@pval[begin[i]:end[i]], se=TRUE)[c(1,2)]
    }
  }
  rs <- ri[ ,1]/ri[ ,2]
  
  size <- as.integer(object@pos[end]-object@pos[begin]+1)
  ord <- order(rs, ri[ ,1], nProbes, size, decreasing=TRUE)
  regions <- data.frame(chr=factor(object@chr[begin])[ord],
                        start=object@pos[begin][ord],
                        end=object@pos[end][ord],
                        size=size[ord], nProbes=nProbes[ord],
                        ri=round(ri[ ,1], 4)[ord], se=round(ri[ ,2], 4)[ord],
                        rs=round(rs, 4)[ord])

  object@regions <- regions
  object@limit <- limit
  object@minLength <- as.integer(minLength)
  object@maxGap <- maxGap
  
  return(object)
}
)


##################################################
## gsri
##################################################
gsri <- function(pval, grenander=FALSE, se=TRUE, custom=FALSE)  {

  cweight <- rep(1, length(pval))
  res <- fitGSRI(pval, NULL, cweight, length(pval),
                 grenander=grenander, se=se, custom=custom)
  res <- c(res, res[1]*res[3])
  names(res) <- c("GSRI", "se", "n", "nReg")

  return(res)
}


##################################################
## threshold
##################################################
setMethod("threshold", "Les",
          function(object, grenander=FALSE, verbose=FALSE, ...)  {
  
  pval <- object@pval
  nProbes <- length(pval)
  erg <- les:::gsri(pval, grenander, se=FALSE)
  nSigProbes <- erg[1]*nProbes
  nSigLower <- floor(nSigProbes)
  cutoff <- c(NA, sort(object@lambda, decreasing=TRUE))[nSigLower+1]
  ## +1: access right probe due to floor !
  
  if(verbose == TRUE)
    print(sprintf("%g %s%g", nSigLower,
                  "significant probes estimated with limit Lambda>=", cutoff))
  
  object@nSigProbes <- nSigProbes
  object@theta <- cutoff
  
  return(object)
}
)


##################################################
## plot
##################################################
setMethod("plot", "Les",
          function(x, y, chr, region=FALSE,
                   xlim, ylim=c(0, 1), error="none",
                   probePch=20, probeCol="black",
                   sigPch=20, sigCol="red3", errorCol="azure4",
                   regionCol=gray, 
                   rug=FALSE, rugSide=1, limit=TRUE,
                   main=NULL, ...)  {

  if(missing(chr))  {
    if(x@nChr == 1)  {
      chr <- levels(x@chr)
      indChr <- rep(TRUE, length(x@pos))
    }
    else
      stop(c("'chr' must be specified. Possible values are: ",
             paste(levels(x@chr), collapse=", ")))
  }
  else  {
    if(length(chr) == 1 && any(x@chr %in% chr))
      indChr <- x@chr %in% chr
    else
      stop("'chr' must have one match")
  }
  pos <- x@pos[indChr]
  lambda <- x@lambda[indChr]
  if(is.logical(limit) && limit == TRUE && length(x@theta) != 0)
    theta <- x@theta
  else
    limit <- FALSE

  if(region == TRUE && length(x@regions) == 0)  {
    region <- FALSE
    warning("'No 'regions' estimated so far.")
  }

  if(missing(xlim))
    xlim <- range(pos)
  ind <- pos >= xlim[1] & pos <= xlim[2]

  plot(pos[ind], lambda[ind], type="n",
       xlim=xlim, ylim=ylim,
       xlab="Probe position", ylab=expression(Lambda), main=main)
  abline(h=c(0,1), col="lightgray")

  val <- pmatch(error, c("se", "ci"), NA)
  if(!is.na(val) && length(x@ci) != 0)  {
    sfrac <- min(min(diff(pos[ind]))/(xlim[2]-xlim[1])/2, 0.01)
    suppressWarnings({
      ss <- ind2log(x@subset, length(x@pos))
      ci <- x@ci[indChr, ]
      gplots::plotCI(pos[ss[indChr]], lambda[ss[indChr]],
                     ui=ci$upper, li=ci$lower,
                     gap=0, pch=".", col=errorCol, add=TRUE,
                     sfrac=sfrac)
    })
  }
  
  points(pos[ind], lambda[ind], type="o", pch=probePch)
  if(limit == TRUE)  {
    sig <- lambda[ind] >= theta
    abline(h=theta, col="gray")
    points(pos[ind][sig], lambda[ind][sig], pch=sigPch, col=sigCol)
  }
  if(rug == TRUE)
    rug(pos[ind], side=rugSide)

  if(region == TRUE && length(x@regions) != 0 && nrow(x@regions) != 0)  {
    regions <- x@regions
    regions <- regions[regions$chr %in% chr, ]
    indRegion <- (regions$start >= xlim[1] & regions$start <= xlim[2]) | (regions$end >= xlim[1] & regions$end <= xlim[2])
    regions <- regions[indRegion, ]
    yr <- c(par()$usr[3], min(ylim)-par()$usr[3])
    if(is.vector(regionCol))
    rect(regions$start, yr[1]+0.2*yr[2], regions$end, yr[1]+0.8*yr[2],
         col=rep(regionCol, length.out=nrow(regions)))    
    else
    rect(regions$start, yr[1]+0.2*yr[2], regions$end, yr[1]+0.8*yr[2],
         col=regionCol(1-regions$ri))
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
                "with window size", object@win))
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


setMethod("export", "Les",
          function(object, file, format="bed", chr, range,
                   description="Lambda", strand=".", group="les",
                   precision=4, ...)  {

  choice <- pmatch(format, c("gff", "bed", "wig"))          
  if(is.na(choice))
    stop("'format' must be 'gff', 'bed' or 'wig'")

  if(choice %in% 1:2)  {
    if(length(object@regions) == 0)
      stop("'regions()' must be run first.")
    regions <- object@regions
    if(nrow(regions) == 0)  {
      warning("No regions to export, no output file written.")
      choice <- 0
    }
    if(missing(chr))
      ind <- rep(TRUE, nrow(regions))
    else
      ind <- regions$chr %in% chr
    chrom <- paste("chr", regions$chr[ind], sep="")
    header <- sprintf("%s%s%s", "track name=LES description=\"",
                      description, "\" useScore=1")
    score <- round(regions$ri[ind]*1e3)
  
    if(choice == 1)  ## gff
      df <- data.frame(chrom=chrom, source="les", feature="region",
                       start=regions$start[ind], end=regions$end[ind],
                       score=score, strand=strand, frame=".", group=group)
    if(choice == 2)  ## bed
      df <- data.frame(chrom=chrom, chromStart=regions$start[ind],
                       chromEnd=regions$end[ind], name=which(ind),
                       score=score)
  }
  if(choice == 3)  {
    if(missing(chr) || length(chr) > 1)
      stop(c("'chr' must have one value for 'wig' export. Possible values are: ",
             paste(levels(object@chr), collapse=", ")))
    header <- c(sprintf("%s%s%s","track name=\"", description,
                        "\" type=wiggle_0 viewLimits=0:1 autoScale=off"),
                sprintf("%s%s %s", "variableStep chrom=", chr,  "span=1"))

    if(missing(range))
      ind <- object@chr %in% chr & !is.na(object@lambda)
    else  {
      ind <- object@chr %in% chr & object@pos >= min(range) &
      object@pos <= max(range) & !is.na(object@lambda)
    }
    
    values <- format(round(object@lambda[ind], precision), scientific=16)
    df <- data.frame(pos=object@pos[ind], value=values)
    indDup <- duplicated(df$pos)
    df <- df[!indDup, ]
    indValid <- !is.na(df$value)
    df <- df[indValid, ]
  }

  if(choice != 0)  {
    write(header, file)
    utils::write.table(df, file, append=TRUE, quote=FALSE, sep="\t",
                       row.names=FALSE, col.names=FALSE)
  }
}
)


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

