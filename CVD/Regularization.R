# TODO: Add comment
# 
# Author: tao.xu
###############################################################################
metabo.asso = scan(what = character())
M22030
M35675
M17945
M32654
M32675
M15122
M32346
M39379
M35126
M35127
M37058
M32319
M16634
M22649
M24074
M32634
M32740
M33132
M33638
M34123
M35193
M35464


###data
tmp = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		scale(S4$ltalteru), scale(S4$ltbmi), as.factor(S4$lcsex), ##model1
		S4$my.diab, ##model 2
		scale(S4$ltsysmm), S4$my.cigreg, S4$my.alkkon, scale(S4$ll_hdla), scale(S4$ll_chola), ##model 3
		scale(S4$lh_crp),
		scale(log(S4[, metabo.asso])), scale(S4[, metabo.ratio.asso])
)
clinical = c("age", "ltbmi", "sex", "diabetes", "smoking", "alkkon","ltsysmm", "ll_hdla", "ll_chola", "lh_crp")
colnames(tmp)[3:12] = clinical#, "total2HDL"
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
		Surv(time, event)~., data = tmp[subset,c("event", "time", metabo.selected3)], 
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


cat("~",paste(metabo.asso, collapse = "+"))

tmp2=tmp[subset,]
model.penal.opt =optL1(
		Surv(time, event), 
		penalized = tmp2[,c(metabo.asso, metabo.ratio.asso)],
		unpenalized = ~ age +ltbmi + sex + diabetes + ltsysmm + ll_hdla + ll_chola + smoking + alkkon +lh_crp,
		data = tmp2,
		fold = 10,	
		minlambda1 = 0.5, maxlambda1 = 10, 
		standardize = T
)
#coxph(Surv(time, event)~., data = tmp[subset, c(1:2, 11:18)])

		
unpenalized = ~ ltalteru +ltbmi + lcsex + lp_diab_who06 + ltsysmm + ll_hdln + ll_choln + strata(ltcigreg) + ltalkkon +lh_crp + total2HDL,

		
metabo.selected2 = scan(what = character())
Arg
PC_ae_C38_0
PC_aa_C32_2
lysoPC_a_C18_1
Trp



metabo.selected3 = scan(what = character())
Arg.Trp
PC_aa_C32_2
lysoPC_a_C17_0


lysoPC_a_C16_0
PC_aa_C28_1
PC_aa_C32_2
PC_aa_C34_2
PC_ae_C40_1
Arg.Trp
Arg.lysoPC_a_C17_0
Arg.lysoPC_a_C18_2 
 

Arg
Trp
lysoPC_a_C17_0
PC_aa_C32_2


Arg
Trp
lysoPC_a_C17_0
PC_aa_C32_2
PC_aa_C34_2
PC_ae_C36_2

Arg
Trp
lysoPC_a_C17_0
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2 
PC_aa_C34_2
PC_aa_C36_3
PC_ae_C38_2
PC_ae_C40_1


C16_1
Arg
Trp
lysoPC_a_C18_1
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2
PC_aa_C34_2
PC_aa_C36_3
PC_ae_C38_0
PC_ae_C38_2
PC_ae_C40_1

Models <- list(
		"Cox.metabolites" = coxph(Surv(time, event) ~ ., data = tmp[subset, c("event", "time", metabo.selected3)]),
		"Cox.clinical" = coxph(Surv(time, event) ~ ., data = tmp[subset,  c("event", "time",clinical)]),
		"Cox.all" = coxph(Surv(time, event) ~ ., data = tmp[subset, c("event", "time", metabo.selected2, clinical)])
)

#### stepwise selection of cox regression
selectCox <- function(formula, data, rule = "aic") {
	require("rms")
	require("prodlim")
	fit <- cph(formula, data, surv = TRUE)
	bwfit <- fastbw(fit, rule = rule)
	if (length(bwfit$names.kept) == 0) {
		newform <- reformulate("1", formula[[2]])
		newfit <- prodlim(newform, data = data)
	} else{
		newform <- reformulate(bwfit$names.kept, formula[[2]])
		newfit <- cph(newform, data, surv = TRUE)
	}
	out <- list(fit = newfit,In = bwfit$names.kept)
	out$call <- match.call()
	class(out) <- "selectCox"
	out
}

selectCox(
		formula = Surv(mi_time, inz_mi) ~  .,
		data = data.frame(tmp[,c(metabo.selected3, clinical)], mi_time = S4$mi_time, inz_mi = S4$inz_mi)[subset, ],
		rule = "p")

+ ltalteru + log(ltdiamm) + log(ltsysmm) + log(ll_hdln) + log(ll_choln) + as.factor(lp_diab_who06) + as.factor(lcsex) + as.factor(ltcigreg)

#sapply(names(Predicts), function(x)  substitute(x))
#(AUC: k), list(k = round(auc[[x]],3))))
#
#c(
#		substitute("all" (AUC: k), list(k = round(auc[[1]],2))),
#		substitute(metabolites (AUC: k), list(k = round(auc[[2]],2))),
#		substitute(clinical (AUC: k), list(k = round(auc[[3]],2)))
#)
#
#pred = prediction(Cox.all$cv.fit, tmp$event[subset])#combined model blue
#plot(performance(pred, "tpr", "fpr"), col = "blue");abline(0,1)
#auc[1] = performance(pred, "auc")@y.values
#pred = prediction(Cox.clinical$cv.fit, tmp$event[subset])#clinical model green
#plot(performance(pred, "tpr", "fpr"), col = "green", add = T)
#auc[2] = performance(pred, "auc")@y.values
#pred = prediction(Cox.metabolites$cv.fit, tmp$event[subset])#metabolites model red
#plot(performance(pred, "tpr", "fpr"), col = "red", add = T)
#auc[3] = performance(pred, "auc")@y.values
#legend(0.6, 0.2, 
#		legend = c(
#				substitute(all (AUC: k), list(k = round(auc[[1]],2))),
#				substitute(metabolites (AUC: k), list(k = round(auc[[2]],2))),
#				substitute(clinical (AUC: k), list(k = round(auc[[3]],2)))
#				), 
#		col = c("blue","red", "green"), lty = 1)
#dev.off()
#auc


#test = function (response, penalized, unpenalized, minlambda1, maxlambda1, 
#		base1, lambda2 = 0, fusedl = FALSE, positive = FALSE, data, 
#		model = c("cox", "logistic", "linear", "poisson"), startbeta, 
#		startgamma, fold, epsilon = 1e-10, maxiter = Inf, standardize = FALSE, 
#		tol = .Machine$double.eps^0.25, trace = TRUE) 
#{
#	prep<-.checkinput(match.call(), parent.frame())
#	#fit <- .modelswitch(prep$model, prep)
#}
#
#test(Surv(time, event)~., unpenalized = ~ ltalteru +ltbmi + lcsex + lp_diab_who06 + ltsysmm + ll_hdln + ll_choln + ltcigreg, data = tmp[subset, c("event", "time", metabo.asso, clinical)],
#		fold = 10,	
#		standardize = T
#)