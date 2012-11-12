# TODO: Add comment
# 
# Author: tao.xu
###############################################################################

#################	change of covariate coefficients while adding metabolites into the model	##################

model = coxph(Surv(mi_time, inz_mi) ~ log(Arg) + log(Trp) + log(lysoPC_a_C17_0) + log(PC_aa_C32_2) +
				ltalteru + factor(lcsex, ordered = F) + ltbmi## model 1
				+ factor(lp_diab_who06==4|lp_diab_who06==5, ordered = F)  ##model 2
				+ log(ltsysmm)+ my.cigreg + my.alkkon  + my.chola + my.hdla ##model 3+ total2HDL
		+ lh_crp  ##model 4 + total2HDL
		#+ltdiamm ## model 5
		#+lh_crp #model 6
		#+waist2hip#model 7
		,subset = which(S4$prev_mi == 0),
		S4)


###################	evaluation by prediction error curve #############
test = sample(subset, 300)
train = setdiff(subset,test)

train = subset

Models <- list(
		"Cox.metabolites" = coxph(Surv(time, event) ~ ., data = tmp[train, c(1:2, 11:19)]),
		"Cox.clinical" = coxph(Surv(time, event) ~ ., data = tmp[train, c(1:10)]),
		"Cox.all" = coxph(Surv(time, event) ~ ., data = tmp[train, -c(17:19)])
)

f = as.formula(paste("Surv(time, event)", paste(colnames(tmp)[-c(1:2)], collapse = "+"), sep = "~"))

PredError <- pec::pec(object = Models, 
		formula = f,
		data = tmp[test,],
		cens.model="marginal",
		splitMethod = "none",
)
plot(PredError, xlim = c(3000, 4000), smooth = T)		

##survival ROC by survivalROC (PS: can only evaluate the ROC of single marker) 
#should not run!
#cutoff = max(tmp$time, na.rm = T)
#nobs = NROW(tmp)
#survival.all = survivalROC.C(
#		Stime = tmp$time[subset],
#		status = tmp$event[subset],
#		marker = tmp[subset, metabo.selected2],
#		predict.time = cutoff,
#		span = 0.25*nobs^(-0.20)
#)



#testing the differences of predicting performances by auc (package: survAUC)
require(survAUC)
Surv.rsp <- Surv(tmp$time[train], tmp$event[train])
Surv.rsp.new <- Surv(tmp$time[test], tmp$event[test])

lpnew = predict(Models$Cox.all, newdata = tmp[test,-c(1:2)])
lp = predict(Models$Cox.all)
times = seq(0,3700, 100)
AUC_sh.all = AUC.sh(Surv.rsp,Surv.rsp.new,lp,lpnew,times)
plot(AUC_sh.all, ylim = c(0.5, 1.0), add = T)
#,ylim = c(0.6, 0.8)
lpnew = predict(Models$Cox.clinical, newdata = tmp[test,-c(1:2)])
lp = predict(Models$Cox.clinical)
times = seq(0,3700, 100)
AUC_sh.clinical = AUC.sh(Surv.rsp,Surv.rsp.new,lp,lpnew,times)
plot(AUC_sh.clinical, add = T, col = "green")

lpnew = predict(Models$Cox.metabolites, newdata = tmp[test,-c(1:2)])
lp = predict(Models$Cox.metabolites)
times = seq(0,3700, 100)
AUC_sh.metabolites = AUC.sh(Surv.rsp,Surv.rsp.new,lp,lpnew,times)
plot(AUC_sh.clinical, add = T, col = "black")

#model comparison by likelihood ratio test
pchisq(-2*(Models$Cox.all$loglik[2] - Models$Cox.clinical$loglik[2]), df=1) 
#model comparison by likelihood ratio test (using anova)
anova(Models$Cox.all, Models$Cox.clinical, Models$Cox.metabolites)

#testing the differences of prediction performance by roc (package: pROC)
set.seed(10)
Predicts <-list(
		#"metabolites.boost" = crossval.cox(x = tmp[subset, metabo.selected], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 1330),
		"metabolites.regularize" = crossval.cox(x = tmp[subset, metabo.selected3], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 1315),
		"clinical" = crossval.cox(x = tmp[subset, clinical], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 1315),
		#"all.boost" = crossval.cox(x = tmp[subset, c(clinical,metabo.selected)], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 1330),	
		"all.regularize" = crossval.cox(x = tmp[subset, c(clinical,metabo.selected3)], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = 1315)
)

auc = NULL
for (i in 1:length(Predicts)){
	#pred = prediction(Predicts[[i]]$cv.fit, tmp$event[subset])
	pred = roc(tmp$event[subset], Predicts[[i]]$cv.fit, ci = T)
	if(i ==1) add = F 
	else add = T
#	plot(performance(pred, "sens", "fpr"), col = i, add = add, lty = i)
#	auc[i] = performance(pred, "auc")@y.values
	plot(pred, col = i, add = add, lty = i)
	auc = rbind(auc, pred$ci[1:3])
	#pred.ci.se = ci.se(pred)
	#plot(pred.ci.se, type = "shape", add = T, col = i)
}
names(auc) = names(Predicts)
#abline(0,1)
legend(0.5, 0.2, 
		legend = paste(names(Predicts), sapply(auc[,2],round,2), sep = " AUC:"), 
		col = c(1:length(Predicts)), lty = c(1:length(Predicts))
)
#pred = prediction(1-survival(model.penal.opt$predictions, 3671), tmp$event[subset])
#plot(performance(pred, "tpr", "fpr"), col = 5, add = T)

pred.1 = roc(tmp$event[subset], Predicts[[3]]$cv.fit, ci = T)
pred.2 = roc(tmp$event[subset], Predicts[[1]]$cv.fit, ci = T)
roc.test(pred.1, pred.2)

#testroc = roc(tmp$event[subset], Predicts[[1]]$cv.fit, ci = T)
#remove(testroc.ci.sp, testroc.ci.threshold, testroc.ci.se, testroc)