estimate <- function(pos, pval, win,
                     weighting=triangWeight, grenander=FALSE, nCores=FALSE)  {

  #checkInputs(pos, pval, win)
  win <- as.integer(win)

  indValid <- !is.na(pval)
  pos <- pos[indValid]
  pval <- round(pval[indValid], 12)

  ord <- sort.list(pos, method="quick", na.last=NA)
  pos <- pos[ord]
  pval <- pval[ord]

  indProbes <- seq(win+1, win+length(pos))
  pos <- c(rep(-Inf, win), pos, rep(Inf, win))
  pval <- c(rep(NA, win), pval, rep(NA, win))

  if(is.numeric(nCores))  {
   cs <- mcsapply(indProbes, calcSingle, pos, pval,
                   win, weighting, grenander, mc.cores=nCores)
  }
  else  {
    cs <- sapply(indProbes, calcSingle, pos, pval,
                 win, weighting, grenander)
  }

  res <- new(Class="Les", pos=cs[1, ], q=cs[2, ],
             se=cs[3, ], nUsedProbes=as.integer(cs[4, ]),
             pval=pval[indProbes], win=win,
             weighting=weighting, grenander=grenander)
  
  return(res)
}


calcSingle <- function(ind0, pos, pval, win,
                       weighting, grenander)  {

  pos0 <- pos[ind0]
  indCut <- seq(ind0-win, ind0+win)
  pos <- pos[indCut]
  pval <- pval[indCut]
  
  distance <- abs(pos - pos0)
  indValid <- distance <= win
  #nValidProbes <- sum(indValid)
  nValidProbes <- length(unique(pval[indValid]))
  #nUniqueProbes <- length(unique(pval[indValid]))

  if(nValidProbes>1)  {
    dis <- distance[indValid]
    weight <- weighting(dis, win)

    res <- fitGSRI(pval[indValid], weight,
                   pos0, nValidProbes, grenander)
  }
  else  {
    res <- c(pos0, NA, NA, nValidProbes)
  }

  return(res)
}


fitGSRI <- function(pval, weight, pos0,
                    nValidProbes, grenander)  {
  
  maxIter <- nValidProbes
  restOld <- 0
  rest <- restOld
  q <- 1

  cdf <- getCDF(pval, weight, grenander)

  for(i in 1:maxIter)  {
    rest <- nValidProbes - ceiling(q*nValidProbes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(c(nValidProbes-1, rest))
    if(is.na(rest) || restOld == rest)
      break
    x <- cdf$sortedPval[rest:nValidProbes] - 1
    y <- cdf$cdf[rest:nValidProbes] - 1
    q <- GSRI:::slopeFast(x, y)
    restOld <- rest
  }
  res <- c(pos0, 1 - min(q, 1), seFast(x, y, q), nValidProbes)

  return(res)
}


getCDF <- function(pval, weight, grenander)  {
  
  ord <- sort.list(pval, method="quick", na.last=NA)
  pval <- pval[ord]
  weight <- weight[ord]

  indUnique <- !duplicated(pval)
  uniquePval <- pval[indUnique]
  uniqueWeight <- weight[indUnique]

  nUnique <- length(uniquePval)
  nProbes <- length(pval)

  cdf <- cumsum(uniqueWeight)
  cdf <- cdf/cdf[nUnique]
  
  if(nProbes != nUnique)
    cdf <- rep.int(cdf, table(pval))
  if(grenander == TRUE)  {
    if(nProbes != nUnique)
      stop("Grenander estimator does not allow duplicates in p-values.")
    cdf <- GSRI:::grenanderInterpol(pval, cdf)
  }
  cdf <- cdf - 0.5/nProbes  ## where to put this ??
  res <- list(sortedPval=pval, cdf=cdf)

  return(res)
}


seFast <- function(x, y, b) {
  
  si     <- sum((y-b*x)^2)/(length(x)-1)
  se     <- as.numeric(sqrt(si/(t(x)%*%x)))
  
  return(se)
}
  

triangWeight <- function(distance, win)  {
  
  weight <- 1 - distance/win  ## not normed
  ## does not matter, normed in cdf
    
  return(weight)
}


rectangWeight <- function(distance, win)  {

  n <- length(distance)
  weight <- rep(1/n, n)

  return(weight)
}


bootSingle <- function(ind0, pos, pval, win, weighting,
                       nBoot, conf, grenander)  {

  pos0 <- pos[ind0]
  indCut <- seq(ind0-win, ind0+win)
  pos <- pos[indCut]
  pval <- pval[indCut]
  
  distance <- abs(pos - pos0)
  indValid <- distance <= win
  #nValidProbes <- sum(indValid)
  nValidProbes <- length(unique(pval[indValid]))

  if(nValidProbes>1)  {
    dis <- distance[indValid]
    weight <- weighting(dis, win)

    bo <- boot(pval[indValid], fitGSRIboot, nBoot,
               cweight=weight, grenander=grenander,
               nValidProbes=nValidProbes)
    ## try not needed?
    ci <- try(boot.ci(bo, conf, type="perc"), silent=FALSE)
    ci <- boot.ci(bo, conf, type="perc")
    if(class(ci)=="try-error")
      res <- c(NA, NA)
    else
      res <- ci$perc[4:5]
  }
  else  {
    res <- c(NA, NA)
  }
  
  return(res)
}


fitGSRIboot <- function(pval, index, cweight,
                        nValidProbes, grenander)  {

  pval <- pval[index]
  weight <- cweight[index]
  maxIter <- nValidProbes
  restOld <- 0
  q <- 1

  cdf <- getCDF(pval, weight, grenander)

  for(i in 1:maxIter)  {
    rest <- nValidProbes - ceiling(q*nValidProbes)
    rest <- max(c(restOld, rest, 1))
    rest <- min(c(nValidProbes-1, rest))
    if(is.na(rest) || restOld == rest)
      break
    x <- cdf$sortedPval[rest:nValidProbes] - 1
    y <- cdf$cdf[rest:nValidProbes] - 1
    if(length(unique(x)) == 1)
      break
    q <- GSRI:::slopeFast(x, y)
    restOld <- rest
  }
  if(is.na(q) || q > 1)  ## good? NA seems not thrown out in boot.ci
    q <- NA
  res <- 1 - min(q, 1)
  
  return(res)
}


mcsapply <- function(X, FUN, ...,
                     simplify=TRUE, USE.NAMES=TRUE, mc.cores)  {
  
  FUN <- match.fun(FUN)
  answer <- mclapply(X, FUN, ..., mc.cores=mc.cores)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
    length)))) == 1L) {
    if (common.len == 1L) 
      unlist(answer, recursive = FALSE)
    else if (common.len > 1L) 
      array(unlist(answer, recursive = FALSE), dim = c(common.len, 
        length(X)), dimnames = if (!(is.null(n1 <- names(answer[[1L]])) & 
        is.null(n2 <- names(answer)))) 
      list(n1, n2))
    else answer
  }
  else answer
}


setClass("Les",
  representation(pos="numeric", q="numeric", se="numeric",
  nUsedProbes="integer", pval="numeric", win="integer",
  weighting="function", grenander="logical", ci="matrix", nBoot="integer", conf="numeric",
  subset="numeric", cutoff="numeric", nSigProbes="numeric",
  regions="matrix", limit="numeric")
)


setGeneric("ci",
           function(res, ...) {standardGeneric("ci")}
           )

setMethod("ci", "Les",
          function(res, subset, nBoot=100, conf=0.95, nCores=FALSE)  {

            if(missing(subset))
              subset <- seq(along=res@pos)
            if(is.logical(subset))
              subset <- which(subset)
            nBoot <- as.integer(nBoot)
            if(res@grenander == TRUE)
              stop("Bootstrap and Grenander estimator do not work together.")

            win <- res@win
            indProbes <- seq(win+1, win+length(pos))[subset]
            pos <- c(rep(-Inf, win), res@pos, rep(Inf, win))
            pval <- c(rep(NA, win), res@pval, rep(NA, win))

            if(is.numeric(nCores))  {
              bs <- mcsapply(indProbes, bootSingle,
                             pos, pval, win, res@weighting,
                             nBoot, conf, res@grenander,
                             mc.cores=nCores)
            }
            else  {
              bs <- sapply(indProbes, bootSingle,
                           pos, pval, win, res@weighting,
                           nBoot, conf, res@grenander)
            }
            
            
            bs <- mcsapply(indProbes, bootSingle,
               pos, pval, win, res@weighting,
               nBoot, conf, res@grenander,
                           mc.cores=multicore:::detectCores()-1)
            
            rownames(bs) <- c("lower", "upper")
            colnames(bs) <- pos[indProbes]

            res@ci <- bs
            res@subset <- as.integer(subset)
            res@nBoot <- nBoot
            res@conf <- conf

            return(res)
}
)


setGeneric("regions",
           function(res, ...) {standardGeneric("regions")}
           )

setMethod("regions", "Les",
          function(res, limit=NULL, verbose=FALSE)  {

            if(is.null(limit))  {
              if(length(res@cutoff) == 0)
                stop("'limit' must be specified")
              else
                limit <- res@cutoff
            }
             
            ind <- c(0, res@q >= limit, 0)
            d <- diff(ind)
            begin <- which(d == 1)
            end <- which(d == -1) - 1
            n <- length(begin)
            if(n != length(end))
              stop("Unequal numbers of region bounderies")
            if(verbose == TRUE)
              print(sprintf("%d %s%g", n,
                            "regions found with q>=", limit))
            regions <- matrix(c(begin, end), 2, n, byrow=TRUE)
            rownames(regions) <- c("begin", "end")
            colnames(regions) <- seq(1, length=n)

            res@regions <- regions
            res@limit <- limit

            return(res)
          }
)


setGeneric("cutoff",
           function(res, ...) {standardGeneric("cutoff")}
           )

setMethod("cutoff", "Les",
          function(res, grenander=FALSE, verbose=FALSE)  {

            pval <- res@pval
            nProbes <- length(pval)
            erg <- fitGSRI(pval, rep(1, nProbes), NA,
                           nProbes, grenander)
            nSigProbes <- erg[2]*nProbes
            nSigLower <- ceiling(nSigProbes)
            cutoff <- c(NA, sort(res@q, decreasing=TRUE))[nSigLower+1]
            ## +1: access right probe due to ceiling !

            if(verbose == TRUE)
              print(sprintf("%g %s%g", nSigLower,
                  "significant probes estimated with limit q>=", cutoff))

            res@nSigProbes <- nSigProbes
            res@cutoff <- cutoff
  
            return(res)
          }
          )


setMethod("plot", "Les",
          function(x, y, region=FALSE, limit=TRUE,
                   xlim, ylim=c(0,1),
                   error="none", semSpread=1.96,
                   patchCol="gray", borderCol="black",
                   markerCol="red", regionCol,
                   regionLty=2, regionLw=3, regionCex=NULL,
                   probeCol="black", probeCex=0.5, probePch=20)  {

            pos <- x@pos
            q <- x@q

            if(error == "ci" && length(x@ci) == 0)  {
              warning("CI not computed so far")
              error <- "none"
            }
    
            if(!missing(y))  {
              par(mar=c(2.1,4.1,4.1,2.1))
              par(fig=c(0,1,0,0.3))
              plot.new()
              plot.window(x=xlim, y=c(0, 2))
              abline(h=0.5, col="gray32")
              indAnno <- y$pos[1, ] >= xlim[1] & y$pos[2, ] <= xlim[2]
              for(i in which(indAnno))  {
                polygon(c(y$pos[ ,i], rev(y$pos[ ,i])),
                        c(0, 0, 1, 1),
                        col=y$col[i])
              }
              text(x=0.5*(anno$pos[2,indAnno]+anno$pos[1,indAnno]),
                   y=0.5,
                   label=anno$name[indAnno],
                   adj=c(0.5,0.5)
                   )
              text(x=xlim[2], y=0.55,
                   label="DNA", adj=c(0.5,0))

              # for next plot
              par(fig=c(0,1,0.2,1), new=TRUE)
            }
            
            if(missing(xlim))
              xlim <- range(pos)
            ind <- (pos >= xlim[1]) & (pos <= xlim[2])  ## needed for plotting?

            if(is.logical(limit) && limit == TRUE)  {
              if(length(res@limit) != 0)
                limit <- res@cutoff
              else  {
                warning("'cutoff' not estimated")
                limit <- FALSE  ## skip warning??
              }
            }
            if(region == TRUE && length(x@regions) != 0)
              regions <- x@regions
            if(!is.logical(region))
              regions <- region

            sig <- q >= limit
            plot(pos, q, type="n", xlim, ylim,
                 xlab="Probe position", ylab=expression(Lambda))
            if(is.numeric(limit))
              abline(h=limit, col="gray")

            #points(pos, q, col=probeCol, cex=probeCex, pch=probePch)

            #sfrac <- min(0.01*range(xlim), min(diff(pos[ind])*0.5))
            sfrac <- min(min(diff(pos[ind]))/(xlim[2]-xlim[1])/2, 0.01)
            suppressWarnings(
            switch(match(error, c("ci", "se", "none")),
                   {ss <- res@subset
                    ci <- res@ci
                    plotCI(pos[ss], q[ss],
                           ui=ci["upper", ], li=ci["lower", ],
                           gap=0, pch=".", col="azure4",
                           add=TRUE, sfrac=sfrac)
                   },
                   {se <- res@se
                    plotCI(pos, q,
                           uiw=se*semSpread, liw=se*semSpread,
                           gap=0, pch=".", col="azure4",
                           add=TRUE, sfrac=sfrac)
                   }
                   )
                             )
                   
            points(pos, q, type="o",
                   pch=probePch, col=probeCol, cex=probeCex)
            if(is.numeric(limit))  {
            points(pos[sig], q[sig],
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


setMethod("[", "Les",
          function(x, i, j, drop)  {
            value <- slot(x, i)
            return(value)
          }
)


setReplaceMethod("[", "Les",
                 function(x, i, j, value)  {
                   slot(x, i) <- value
                   validObject(x)
                   return(x)
                 }
)


setMethod("show", "Les",
          function(object)  {
            pos <- object@pos
            q <- object@q
            pr <- range(pos)
            qr <- range(q)
            cat("Object of class 'Les'")
            cat(sprintf("%d %s (%d-%d)", length(pos), "probes", pr[1], pr[2]))
            cat(sprintf("%s %g-%g", "Range Lambda:", qr[1], qr[2]))
            if(length(object@ci) != 0)
              cat("confidence intervals computed", quote=FALSE)
            if(length(object@regions) != 0)
              cat("regions detected", quote=FALSE)
            if(length(object@cutoff) != 0)
              cat("cutoff estimated", quote=FALSE)
          }
          )


setMethod("summary", "Les",
          function(object, ...)  {
          slotNames(object)
        }
)




