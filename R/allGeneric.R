##################################################
## generics for class 'Les'
##################################################

setGeneric("ci",
           function(object, subset, nBoot=100, conf=0.95,
                    nCores=NULL, ...)
           {standardGeneric("ci")}
           )

setGeneric("estimate",
           function(object, win, weighting=triangWeight,
                    grenander=TRUE, se=FALSE, minProbes=3,
                    method="la", nCores=NULL, verbose=FALSE, ...)
           {standardGeneric("estimate")}
           )

setGeneric("regions",
           function(object, limit=NULL, minLength=2, maxGap=Inf,
                    verbose=FALSE, ...)
           {standardGeneric("regions")}
           )

setGeneric("threshold",
           function(object, grenander=FALSE, verbose=FALSE, ...)
           {standardGeneric("threshold")}
           )

setGeneric("export",
           function(object, file, format="bed", chr, range,
                    description="Lambda", strand=".", group="les",
                    precision=4, ...)
           {standardGeneric("export")}
           )

setGeneric("chi2",
           function(object, winSize, regions, offset, fdr="lfdr",
                    method, scaling=les:::scaleNorm, nCores=NULL,
                    verbose=FALSE, ...)
           {standardGeneric("chi2")}
           )
