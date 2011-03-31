###################################################
### chunk number 1: 
###################################################
#line 77 "les.Rnw"
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
#line 122 "les.Rnw"
library(limma)
design <- cbind(offset=1, diff=condition)
fit <- lmFit(exprs, design)
fit <- eBayes(fit)
pval <- fit$p.value[, "diff"]


###################################################
### chunk number 4: plotProbeLevelStatistics
###################################################
#line 130 "les.Rnw"
plot(pos, pval, pch=20, xlab="Probe position", ylab=expression(p))
abline(v=region)


###################################################
### chunk number 5: constructLes
###################################################
#line 169 "les.Rnw"
res <- Les(pos, pval)


###################################################
### chunk number 6: estimateLes
###################################################
#line 187 "les.Rnw"
res <- estimate(res, win=200)


###################################################
### chunk number 7: showPlotLes
###################################################
#line 195 "les.Rnw"
res
plot(res)
abline(v=region)


###################################################
### chunk number 8: estimateLes2
###################################################
#line 210 "les.Rnw"
res <- estimate(res, win=200, weighting=rectangWeight)


###################################################
### chunk number 9: showPlotLes2
###################################################
#line 214 "les.Rnw"
res
plot(res)
abline(v=region)


###################################################
### chunk number 10: threshold
###################################################
#line 241 "les.Rnw"
res <- threshold(res, grenander=TRUE, verbose=TRUE)


###################################################
### chunk number 11: regions
###################################################
#line 256 "les.Rnw"
res <- regions(res)
res
res["regions"]


###################################################
### chunk number 12: plotRegions
###################################################
#line 262 "les.Rnw"
plot(res, region=TRUE)
abline(v=region)


###################################################
### chunk number 13: ci
###################################################
#line 283 "les.Rnw"
subset <- pos >= 5232300 & pos <= 5233200
res <- ci(res, subset, nBoot=50)


###################################################
### chunk number 14: plotCi
###################################################
#line 288 "les.Rnw"
plot(res, error="ci")


###################################################
### chunk number 15: plotOptions
###################################################
#line 309 "les.Rnw"
plot(res, error="ci", region=TRUE, rug=TRUE,
xlim=c(5232000, 5233000), plot=list(main="LES for simulated data"), sigArgs=list(pch="*"), limitArgs=list(lty=2, lwd=3), regionArgs=list(col="black", density=20))


###################################################
### chunk number 16: customWeightingFunction
###################################################
#line 340 "les.Rnw"
weightFoo <- function(distance, win) {
weight <- 1 - distance/win
return(weight)
}
resFoo <- estimate(res, 200, weighting=weightFoo)


###################################################
### chunk number 17: sessionInfo
###################################################
#line 357 "les.Rnw"
toLatex(sessionInfo(), locale=FALSE)


