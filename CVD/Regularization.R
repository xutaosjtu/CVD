# TODO: Add comment
# 
# Author: tao.xu
###############################################################################

###data
tmp = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		S4$ltalteru, S4$ltbmi, as.factor(S4$lcsex), ##model1
		(S4$lp_diab_who06==4|S4$lp_diab_who06==5), ##model 2
		log(S4$ltsysmm),log(S4$ll_hdln), log(S4$ll_choln), as.factor(S4$ltcigreg),##model 3
		##log(S4$total2HDL), ##model 4
		##log(S4$ltdiamm), ##model 5
		log(as.matrix(S4[, metabo.asso]))
)
colnames(tmp)[3:10] = clinical#, "total2HDL"
na.index = unique(unlist(apply(tmp, 2, function(x) which(is.na(x)))))
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi)), na.index)

##regularization based selection
require(penalized)
model.penal = penalized(
		Surv(time, event)~., data = tmp[subset,c("event", "time", metabo.asso)], 
		standardize = T, ##centralize the covariates 
		#steps = "Park", trace = F, ## starting from the largest value of lambda1 untile the specificed lambda1 is reached, the function will return a list of penfit objects, and the chang of coefficients could be visualized by plotpath 
		#positive = T, ## positive restriction to all regression coefficients
		lambda1 = 2, lambda2= 0 #penalization parameter 
)
plotpath(model.penal, log = "x")
#plot(coefficients(model.penal,"all"),main = "lasso", col="red",xlab = "probes",ylab = "coefficients",type="l")##plot the coefficients

require(globaltest)##pretesting by global test
gt(Surv(tmp$time[subset], tmp$event[subset]), tmp[subset, c(11:19)])

model.penal = cvl(
		Surv(time, event)~., data = tmp[subset,c("event", "time", metabo.asso)], 
		standardize = T, ##centralize the covariates 
		#steps = "Park", trace = F, ## starting from the largest value of lambda1 untile the specificed lambda1 is reached, the function will return a list of penfit objects, and the chang of coefficients could be visualized by plotpath 
		#positive = T, ## positive restriction to all regression coefficients
		fold = model.penal$fold, ## k-fold cross-validation
		lambda1 = 2, lambda2= 0#penalization parameter
)

model.penal = profL1(
		Surv(time, event)~., data = tmp[subset,c("event", "time", metabo.asso)],
		fold = 10,
		plot = T,
		save.predictions = T,
		minl = 0.01, maxl = 5
)
plot(model.penal$lambda, model.penal$cvl, type = "l", log = "x")
plotpath(model.penal$fullfit, log = "x")
pred = prediction((1-survival(model.penal.opt$predictions,2833)), tmp$event[subset])#ROC evaluation at each time point, after profiling of the model at different lamda
plot(performance(pred, "tpr", "fpr"), add= T, col = "red");abline(0,1)
performance(pred, "auc")

model.penal.opt =optL1(
		Surv(time, event)~., data = tmp[subset, c("event", "time", metabo.asso)],
		fold = model.penal$fold,
		standardize = T
)
coxph(Surv(time, event)~., data = tmp[subset, c(1:2, 11:18)])


metabo.selected2 = scan(what = character())
Arg
PC_ae_C38_0
PC_aa_C32_2
lysoPC_a_C18_1
Trp


Models <- list(
		"Cox.metabolites" = coxph(Surv(time, event) ~ ., data = tmp[subset, c("event", "time", metabo.selected2)]),
		"Cox.clinical" = coxph(Surv(time, event) ~ ., data = tmp[subset,  c("event", "time",clinical)]),
		"Cox.all" = coxph(Surv(time, event) ~ ., data = tmp[subset, c("event", "time", metabo.selected2, clinical)])
)

theta.fit <- function(x, y, ...) {
	d = data.frame(y, x)
	#print(colnames(d))
	coxph(y ~  ., data=d)
}
theta.predict <- function(fit, x){
	#if(is.null(dim(x))) x=t(x)
	#dim(x)
	#print(colnames(x))
	value=predict(fit,newdata= as.data.frame(x) , type="risk")
	return(value)
}

Cox.metabolites = crossval.cox(x = tmp[subset, metabo.selected2], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 10)

x = tmp[subset, c(clinical, metabo.asso)]
y= Surv(tmp$time[subset], tmp$event[subset])
n = dim(x)[1]
ngroup = 10
leave.out <- trunc(n/ngroup)
o <- sample(1:n)
groups <- vector("list", ngroup)
for (j in 1:(ngroup - 1)) {
	jj <- (1 + (j - 1) * leave.out)
	groups[[j]] <- (o[jj:(jj + leave.out - 1)])
}
groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
u = NULL
cv.fit = rep(NA, n)
for (j in 1:ngroup) {
	u <- theta.fit(x[-groups[[j]], ], y[-groups[[j]]])
	cv.fit[groups[[j]]] <- theta.predict(u, x[groups[[j]],])
}
Cox.all2 = list(cv.fit = cv.fit, ngroup = ngroup, leave.out = leave.out, groups = groups)


#Cox.all = crossval.cox(x = tmp[subset, c(clinical, metabo.selected2)], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 10)

auc = vector("numeric", 3)
pred = prediction(Cox.all$cv.fit, tmp$event[subset])#combined model blue
plot(performance(pred, "tpr", "fpr"), col = "blue");abline(0,1)
auc[1] = performance(pred, "auc")@y.values
pred = prediction(Cox.clinical$cv.fit, tmp$event[subset])#clinical model green
plot(performance(pred, "tpr", "fpr"), col = "green", add = T)
auc[2] = performance(pred, "auc")@y.values
pred = prediction(Cox.metabolites$cv.fit, tmp$event[subset])#metabolites model red
plot(performance(pred, "tpr", "fpr"), col = "red", add = T)
auc[3] = performance(pred, "auc")@y.values
legend(0.6, 0.2, 
		legend = c(
				substitute(all (AUC: k), list(k = round(auc[[1]],2))),
				substitute(metabolites (AUC: k), list(k = round(auc[[2]],2))),
				substitute(clinical (AUC: k), list(k = round(auc[[3]],2)))
				), 
		col = c("blue","red", "green"), lty = 1)
dev.off()
auc
