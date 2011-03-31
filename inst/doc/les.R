###################################################
### chunk number 1: 
###################################################
#line 47 "les.Rnw"
set.seed(1)
options(width=65, SweaveHooks=list(fig=function() par(mar=c(5.1, 5.1, 2.1, 1.1))))


###################################################
### chunk number 2: loadData
###################################################
#line 101 "les.Rnw"
library(les)
data(spikeInData)
head(exprs)
dim(exprs)
pos <- as.integer(rownames(exprs))
condition <- as.integer(colnames(exprs))
reference
region <- as.vector(reference[ ,c("start", "end")])


###################################################
### chunk number 3: estimateProbeLevelStatistics
###################################################
#line 118 "les.Rnw"
library(limma)
design <- cbind(offset=1, diff=condition)
fit <- lmFit(exprs, design)
fit <- eBayes(fit)
pval <- fit$p.value[, "diff"]


###################################################
### chunk number 4: plotProbeLevelStatistics
###################################################
#line 126 "les.Rnw"
plot(pos, pval, pch=20, xlab="Probe position", ylab=expression(p))
abline(v=region)


###################################################
### chunk number 5: constructLes
###################################################
#line 165 "les.Rnw"
res <- Les(pos, pval)


###################################################
### chunk number 6: estimateLes
###################################################
#line 177 "les.Rnw"
res <- estimate(res, win=200)


###################################################
### chunk number 7: showPlotLes
###################################################
#line 186 "les.Rnw"
res
summary(res)
plot(res)
abline(v=region)


###################################################
### chunk number 8: showPlotLes2
###################################################
#line 205 "les.Rnw"
res2 <- estimate(res, win=200, weighting=rectangWeight)
res2
plot(res2)
abline(v=region)


###################################################
### chunk number 9: threshold
###################################################
#line 226 "les.Rnw"
res2 <- threshold(res2, grenander=TRUE, verbose=TRUE)


###################################################
### chunk number 10: regions
###################################################
#line 239 "les.Rnw"
res2 <- regions(res2, verbose=TRUE)
res2
res2["regions"]


###################################################
### chunk number 11: plotRegions
###################################################
#line 245 "les.Rnw"
plot(res2, region=TRUE)
abline(v=region)


###################################################
### chunk number 12: plotCi
###################################################
#line 265 "les.Rnw"
subset <- pos >= 5232400 & pos <= 5233100
res2 <- ci(res2, subset, nBoot=50, alpha=0.1)
plot(res2, error="ci", region=TRUE)


###################################################
### chunk number 13: plotOptions
###################################################
#line 289 "les.Rnw"
plot(res2, error="ci", region=TRUE, rug=TRUE, xlim=c(5232000, 5233000), sigArgs=list(col="firebrick4"), plotArgs=list(main="LES results", yaxp=c(0, 1, 2)), limitArgs=list(lty=2, lwd=3), regionArgs=list(col="black", density=20), probeArgs=list(col="dodgerblue4", type="p"))


###################################################
### chunk number 14: export
###################################################
#line 309 "les.Rnw"
bedFile <- paste(tempfile(), "bed", sep=".")
gffFile <- paste(tempfile(), "gff", sep=".")
wigFile <- paste(tempfile(), "wig", sep=".")
export(res2, bedFile)
export(res2, gffFile, format="gff")
export(res2, wigFile, format="wig")


###################################################
### chunk number 15: customWeightingFunction
###################################################
#line 331 "les.Rnw"
weightFoo <- function(distance, win) {
weight <- 1 - distance/win
return(weight)
}
resFoo <- estimate(res, 200, weighting=weightFoo)


###################################################
### chunk number 16: chi2 eval=FALSE
###################################################
## #line 339 "les.Rnw"
## regions <- res["regions"]
## winsize <- seq(100, 300, by=20)
## res2 <- chi2(res2, winsize, regions, offset=2500)
## plot(winsize, x["chi2"], type="b")


###################################################
### chunk number 17: sessionInfo
###################################################
#line 355 "les.Rnw"
toLatex(sessionInfo(), locale=FALSE)


