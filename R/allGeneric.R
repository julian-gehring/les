##################################################
## generics for class 'Les'
##################################################

setGeneric("ci",
           function(object, subset, nBoot=100, conf=0.95,
                    nCores=NULL)
           {standardGeneric("ci")}
           )

setGeneric("estimate",
           function(object, win, weighting=triangWeight,
                    grenander=FALSE, se=TRUE, nCores=NULL)
           {standardGeneric("estimate")}
           )

setGeneric("regions",
           function(object, limit=NULL, minLength, maxGap,
                    verbose=FALSE)
           {standardGeneric("regions")}
           )

setGeneric("threshold",
           function(object, grenander=FALSE, verbose=FALSE)
           {standardGeneric("threshold")}
           )

setGeneric("exportLambda",
           function(object, file, chr, range,
                    description="Lambda", precision=4)
           {standardGeneric("exportLambda")}
           )
