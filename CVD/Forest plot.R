### load metafor package
library(metafor)

### load BCG vaccine dataset
data(dat.bcg)

### calculate log relative risks and corresponding sampling variances
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)

### fit random-effects models to some subsets
res.r1 <- rma(yi, vi, data=dat[1:3,], method = "FE")
res.r2 <- rma(yi, vi, data=dat[1:3,], method = "REML")
res.r3 <- rma(yi, vi, data=dat[1:3,], method = "DL")

### collect model estimates and corresponding variances
estimates <- c(res.r1$yi, res.r2$yi, res.r3$yi)
variances <- c(res.r1$vi, res.r2$vi, res.r3$vi)

### create vector with labels
labels <- rep(c("Study 1", "Study 2", "Study 3"), 3)

##########################
##	First Solution		##
##########################
### forest plot
forest(estimates, variances, slab = labels, xlab="Relative Risk (log scale)", rows = c(3:5, 10:12, 17:19), ylim=c(-1, 24), cex = 0.75)
## you can add polygons to the meta analysis results
addpoly(coef(res.r1), vcov(res.r1), rows=1, col="black", mlab="FE Method", cex = 0.75)
addpoly(coef(res.r2), vcov(res.r2), rows=8, col="black", mlab="REML Method", cex = 0.75)
addpoly(coef(res.r3), vcov(res.r3), rows=15, col="black", mlab="DL Method", cex = 0.75)

##########################
##	Second Solution		##
##########################
### forest plot
par(mfrow=c(1,3))
forest(res.r1, annotate=F, cex = 0.75)
forest(res.r2, annotate=F, slab = rep(NA, 3), mlab = NA, cex = 0.75)
forest(res.r3, annotate=F, slab = rep(NA, 3), mlab = NA, cex = 0.75)