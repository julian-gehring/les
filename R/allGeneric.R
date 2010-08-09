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
                    grenander=FALSE, se=FALSE, minProbes=3,
                    method="la", nCores=NULL, ...)
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
