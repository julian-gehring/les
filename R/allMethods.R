######################################################################
## methods for les package
######################################################################

######################################################################
## estimate
######################################################################
setMethod("estimate", "Les",
          function(object, win, weighting=triangWeight, grenander=TRUE, 
                   se=FALSE, minProbes=3, method="la", nCores=NULL,
                   verbose=FALSE, ...)  {   

  ## check input
  les:::checkState(object@state, "Les", "'Les'")
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

  ## for each chr
  for(c in 1:object@nChr)  {
    if(verbose == TRUE)
      print(sprintf("%s '%s' (%d/%d)",
                    "Chromosome", chrLevel[c], c, object@nChr))
    indChr <- object@chr == chrLevel[c]
    indProbes <- seq(win+1, win+length(object@pos[indChr]))
    pos <- c(rep(-Inf, win), object@pos[indChr], rep(Inf, win))
    pval <- c(rep(NA, win), object@pval[indChr], rep(NA, win))

    ## apply
    cs <- les:::mcsapply(indProbes, les:::calcSingle, pos, pval,
                         win, weighting, grenander, se, nBoot=FALSE,
                         minProbes=minProbes, custom=custom, mc.cores=nCores)

    ## extract result
    object@lambda[indChr] <- cs[1, ]
    object@nProbes[indChr] <- as.integer(cs[3, ])
    if(se == TRUE)
      object@se[indChr] <- cs[2, ]
    if(grenander == TRUE)
      object@lambda0[indChr] <- cs[4, ]
  }

  object@ci <- data.frame()
  object@theta <- numeric()
  object@regions <- data.frame()
  object@nBoot <- integer()
  object@conf <- numeric()
  object@subset <- integer()  

  object@win <- win
  object@weighting <- weighting
  object@grenander <- grenander
  object@minProbes <- minProbes
  object@method <- method
  object@state <- les:::setState(object@state, "estimate")

  return(object)
}
)


######################################################################
## ci
######################################################################
setMethod("ci", "Les",
          function(object, subset, nBoot=100, conf=0.95,
                   nCores=NULL, ...)  {

  les:::checkState(object@state, "estimate")
  if(missing(subset))
    subset <- rep(TRUE, length(object@pos))
  if(!is.logical(subset) & is.vector(subset))
    subset <- ind2log(subset, length(object@pos))
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
    bs <- les:::mcsapply(indProbes, les:::calcSingle, pos, pval,
                         win, object@weighting, object@grenander, se=FALSE,
                         nBoot, conf, minProbes=object@minProbes, custom=custom,
                         mc.cores=nCores)
    bc <- cbind(bc, bs)
  }
  
  object@ci <- data.frame(lower=bc[1, ], upper=bc[2, ])
  object@subset <- log2ind(subset)
  object@nBoot <- nBoot
  object@conf <- conf
  object@state <- les:::setState(object@state, "ci")
  
  return(object)
}
)


######################################################################
## regions
######################################################################
setMethod("regions", "Les",
          function(object, limit=NULL, minLength=2,
                   maxGap=Inf, verbose=FALSE, ...)  {

  les:::checkState(object@state, "estimate")
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
      ri[i, ] <- les:::gsri(object@pval[begin[i]:end[i]], se=TRUE)[c(1,2)]
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
  object@state <- les:::setState(object@state, "regions")
  
  return(object)
}
)


######################################################################
## threshold
######################################################################
setMethod("threshold", "Les",
          function(object, grenander=FALSE, verbose=FALSE, ...)  {

  les:::checkState(object@state, "Les")
  pval <- object@pval
  nProbes <- length(pval)
  erg <- les:::gsri(pval, grenander, se=FALSE)
  nSigProbes <- erg[1]*nProbes
  nSigLower <- floor(nSigProbes)
  names(nSigLower) <- NULL
  theta <- c(NA, sort(object@lambda, decreasing=TRUE))[nSigLower+1]
  
  if(verbose == TRUE)  {
    if(length(object@lambda) != 0)
      print(sprintf("%g %s%g", nSigLower,
                     "significant probes estimated with limit Lambda>=", theta))
    else
      print(sprintf("%g %s", nSigLower, "significant probes estimated"))
  }
  
  object@nSigProbes <- nSigLower
  object@theta <- theta
  object@state <- les:::setState(object@state, "threshold")
  
  return(object)
}
)


######################################################################
## plot
######################################################################
setMethod("plot", "Les",
          function(x, y, chr, region=FALSE,
                   xlim, ylim=c(0, 1), error="none",
                   probePch=20, probeCol="black", lty=1, 
                   sigPch=20, sigCol="red3", errorCol="azure4",
                   regionCol=gray, 
                   rug=FALSE, rugSide=1, limit=TRUE,
                   main=NULL, ...)  {

  les:::checkState(x@state, "estimate")         
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

  graphics::plot(pos[ind], lambda[ind], type="n",
                 xlim=xlim, ylim=ylim,
                 xlab="Probe position", ylab=expression(Lambda), main=main)
  graphics::abline(h=c(0,1), col="lightgray")

  val <- pmatch(error, c("se", "ci"), NA)
  if(!is.na(val) && length(x@ci) != 0)  {
    sfrac <- min(min(diff(pos[ind]))/(xlim[2]-xlim[1])/2, 0.01)
    suppressWarnings({
      ss <- ind2log(x@subset, length(x@pos))
      ci <- x@ci[indChr, ]
      gplots::plotCI(pos[ss[indChr]], lambda[ss[indChr]],
                     ui=ci$upper, li=ci$lower, gap=0,
                     pch=NA, barcol=errorCol, col=probeCol,
                     add=TRUE, sfrac=sfrac)
    })
  }
  
  graphics::points(pos[ind], lambda[ind], type="o", pch=probePch,
                   col=probeCol, lty=lty)
  if(limit == TRUE)  {
    sig <- lambda[ind] >= theta
    graphics::abline(h=theta, col="gray")
    graphics::points(pos[ind][sig], lambda[ind][sig], pch=sigPch, col=sigCol)
  }
  if(rug == TRUE)
    graphics::rug(pos[ind], side=rugSide)

  if(region == TRUE && length(x@regions) != 0 && nrow(x@regions) != 0)  {
    regions <- x@regions
    regions <- regions[regions$chr %in% chr, ]
    indRegion <- (regions$start >= xlim[1] & regions$start <= xlim[2]) | (regions$end >= xlim[1] & regions$end <= xlim[2])
    regions <- regions[indRegion, ]
    yr <- c(par()$usr[3], min(ylim)-par()$usr[3])
    if(is.vector(regionCol))
    graphics::rect(regions$start, yr[1]+0.2*yr[2], regions$end, yr[1]+0.8*yr[2],
                   col=rep(regionCol, length.out=nrow(regions)))    
    else
    graphics::rect(regions$start, yr[1]+0.2*yr[2], regions$end, yr[1]+0.8*yr[2],
         col=regionCol(1-regions$ri))
  }
}
)


######################################################################
## [
######################################################################
setMethod("[", "Les",
          function(x, i, j, drop)  {
            value <- slot(x, i)
            return(value)
          }
)


######################################################################
## <-
######################################################################
setReplaceMethod("[", "Les",
                 function(x, i, j, value)  {
                   slot(x, i) <- value
                   validObject(x)
                   return(x)
                 }
)


######################################################################
## show
######################################################################
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
  if(length(object@nSigProbes) != 0)  {
    if(!is.na(object@theta) && length(object@lambda) != 0)
      cat(sprintf("* %d %s %g\n", ceiling(object@nSigProbes),
                  "regulated probes estimated for lambda >=",
                  object@theta))
    else
      cat(sprintf("* %d %s\n", ceiling(object@nSigProbes),
                  "regulated probes estimated"))
  }
  if(length(object@regions) != 0)
    cat(sprintf("* %d %s\n", nrow(object@regions), "regions detected"))
}
)


######################################################################
## summary
######################################################################
setMethod("summary", "Les",
          function(object, ...)  {

  show(object, ...)
  cat("* Available slots:\n")
  sn <- slotNames(object)
  ind <- sapply(sn, function(sn) length(object[sn])) > 0
  print(sn[ind])
}
)


######################################################################
## export
######################################################################
setMethod("export", "Les",
          function(object, file, format="bed", chr, range,
                   description="Lambda", strand=".", group="les",
                   precision=4, ...)  {

  les:::checkState(object@state, "estimate")
  choice <- pmatch(format, c("gff", "bed", "wig"))          
  if(is.na(choice))
    stop("'format' must be 'gff', 'bed' or 'wig'")

  if(choice == 3)  {
    if(missing(chr))  {
      if(object@nChr == 1)  {
        chr <- levels(object@chr)
      }
      else  {
        stop(c("'chr' must be specified. Possible values are: ",
               paste(levels(object@chr), collapse=", ")))
      }
    }
    else  {
      if(length(chr) != 1 || !any(object@chr %in% chr))
        stop(c("'chr' must have one match. Possible values are: ",
               paste(levels(object@chr), collapse=", ")))
    }
  }
  
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
    else  {
      if(!any(object@chr %in% chr))
        stop(c("'chr' must have one match. Possible values are: ",
               paste(levels(object@chr), collapse=", ")))
      ind <- regions$chr %in% chr
    }
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


######################################################################
## chi2
######################################################################
setMethod("chi2", "Les",
          function(object, winSize, regions, offset,
                   fdr="lfdr", method, scaling=les:::scaleNorm,
                   nCores=NULL, verbose=FALSE, ...)  {

  ## check input arguments
  les:::checkState(object@state, "Les")
  if(missing(regions))  {
    if(length(object@regions) != 0)  {
      if(missing(offset))
        warning("A 'offset' should be defined if working with estimated regions.")
      regions <- object@regions
    }
    else  {
      stop("'regions' must be specified.")
    }
  }
  if(!is.data.frame(regions) || is.null(regions$start) || is.null(regions$end))  {
    stop("'regions' must be a data frame with variables 'start' and 'end'.")
  }
  if(!missing(offset))  {
    regions$start <- regions$start - offset
    regions$end <- regions$end + offset
  }

  winSize <- as.integer(winSize)
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

  ## construct data object
  if(missing(method))
    method <- object@method

  chi2 <- les:::mcsapply(1:nReg, optimalSingleRegion, regions, object, winSize, fdr, method, scaling, nCores, verbose=verbose)

  rownames(chi2) <- winSize
  colnames(chi2) <- rownames(regions)

  ## unset weighting in options
  options(weighting=oldOpt)

  object@winSize <- winSize
  object@chi2 <- chi2
  object@state <- les:::setState(object@state, "chi2")

  return(object)
}
)
