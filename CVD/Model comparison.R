# TODO: Add comment
# 
# Author: tao.xu
###############################################################################
S4$my.sysmm.untreat = log(S4$ltsysmm)
S4$my.sysmm.untreat[which((S4$ltmbbl ==1 | S4$ltmace ==1 | S4$ltmata ==1))] = 0
S4$my.sysmm.treat = log(S4$ltsysmm)
S4$my.sysmm.treat[which(!(S4$ltmbbl ==1 | S4$ltmace ==1 | S4$ltmata ==1))] = 0
data = S4
data[, metabo.selected3] =  log(S4[, metabo.selected3])

#################	Framingham score	##########################
framingham<-function(x){
	#print(x["lcsex"])
	if(x["lcsex"] == 1){ beta = c(3.06117, 1.12370, -0.93263, 1.99881, 1.93303, 0.65451, 0.57367, 23.9802); S0 = 0.88936}
	else {beta = c(2.32888, 1.20904, -0.70833, 2.76157, 2.82263, 0.52873, 0.69154, 26.1931); S0 = 0.95012}
	X = rep(0, 6)
	X[1] = beta[1] * log(x$ltalteru)#age
	X[2] = beta[2] * log(x$ll_chola)#total cholesterol
	X[3] = beta[3] * log(x$ll_hdla)#HDL
	if(is.na(x$ltmbbl) & is.na(x$ltmace) & is.na(x$ltmata)) { X[4]=NA }
	else if(!(x$ltmbbl ==1 | x$ltmace ==1 | x$ltmata ==1)) X[4]= beta[4]*log(x$ltsysmm) #systolic blood pressure if un-treated
	else  X[4] = beta[5]*log(x$ltsysmm) #systolic blood pressure if treated
	X[5] = beta[6]* as.numeric(x$ltcigreg==1|x$ltcigreg==2)#smoking
	X[6] = beta[7] * (as.numeric(x$my.diab)-1)
	risk = 1-S0^exp(sum(X) - beta[8]);
	#risk = exp(sum(X) - beta[8])
	return(risk)
}

for(i in 1:dim(data)[1]){
	data$framingham[i] = framingham(data[i,]) 
}


prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)
model = coxph(Surv(mi_time, inz_mi) ~ Arg + Trp + lysoPC_a_C17_0 + PC_aa_C32_2 +  
				log(ltalteru) + #age
				log(ll_chola) + #total cholesterol
				log(ll_hdla) + #HDL cholesterol
				my.sysmm.untreat + 	#systolic blood pressure if untreated
				my.sysmm.treat +	#systolic blood pressure if treated
				as.numeric(ltcigreg==1|ltcigreg==2) +
				(as.numeric(my.diab)-1),
		subset = which(S4$prev_mi == 0&S4$lcsex ==1),
		data)
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
model = coxph(Surv(mi_time, inz_mi) ~ Arg + Trp + lysoPC_a_C17_0 + PC_aa_C32_2 +  
				log(ltalteru) + #age
				log(ll_chola) + #total cholesterol
				log(ll_hdla) + #HDL cholesterol
				my.sysmm.untreat + 	#systolic blood pressure if untreated
				my.sysmm.treat +	#systolic blood pressure if treated
				as.numeric(ltcigreg==1|ltcigreg==2) +
				(as.numeric(my.diab)-1),
		subset = which(S4$prev_mi == 0&S4$lcsex ==2),
		data)
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")

subset =  which(S4$prev_mi == 0&S4$lcsex ==2)
pred.3 = roc (data$inz_mi[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)

pred.framingham = roc (data$inz_mi[which(S4$prev_mi == 0)], data$framingham[which(S4$prev_mi == 0)], ci = T)

base = data.frame(ltalteru=0, ll_chola =0, ll_hdla = 0, my.sysmm.untreat = 0, my.sysmm.treat = 0, ltcigreg = 0, my.diab = 0)

#################	change of covariate coefficients while adding metabolites into the model	##################
#Arg = scale(S4$Arg); 
#Trp = scale(S4$Trp);
#lysoPC_a_C17_0 = scale(log(S4$lysoPC_a_C17_0))
#PC_aa_C32_2 = scale(log(S4$PC_aa_C32_2))

prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)
model = coxph(Surv(mi_time, inz_mi) ~ Arg + Trp + lysoPC_a_C17_0 + PC_aa_C32_2 +  
				ltalteru + ltbmi +## model 1
				my.diab +  ##model 2
				ltsysmm+ my.cigreg  + my.alkkon + ll_chola + ll_hdla ##model 3+ my.chola + my.hdla + ll_chola + total2HDL +ll_hdla
				+ lh_crp  ##model 4 + my.physical
		,subset = which(S4$prev_mi == 0&S4$lcsex ==1),
		data)
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
model = coxph(Surv(mi_time, inz_mi) ~ Arg + Trp + lysoPC_a_C17_0 + PC_aa_C32_2 +  
				ltalteru + ltbmi +## model 1
				my.diab +  ##model 2
				ltsysmm+ my.cigreg  + my.alkkon + ll_chola + ll_hdla ##model 3+ my.chola + my.hdla + ll_chola + total2HDL +ll_hdla
				+ lh_crp  ##model 4 + my.physical
		,subset = which(S4$prev_mi == 0&S4$lcsex ==2),
		data)
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")

extractAIC(model)


model = coxph(Surv(mi_time, inz_mi) ~ Arg + Trp + lysoPC_a_C17_0 + PC_aa_C32_2 +  
				ltalteru + as.factor(lcsex) + ltbmi +## model 1
				my.diab +  ##model 2
				ltsysmm+ my.cigreg  + my.alkkon + ll_chola + ll_hdla ##model 3+ my.chola + my.hdla + ll_chola + total2HDL +ll_hdla
				+ lh_crp  ##model 4 + my.physical
		,subset = which(S4$prev_mi == 0),
		data)
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")

tmp.pred = crossval.cox(x = tmp[subset, c(metabo.selected3, clinical)], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = length(subset))



plotCI(1:length(coef(model))+0.5, y = coef(model), liw = coef(model)-confint(model)[,1], uiw = abs(coef(model)-confint(model)[,2]), col = "green")

abline(h = 0)
plotCI(1:length(coef(model)[5:14]), y = coef(model)[5:14], liw = coef(model)[5:14]-confint(model)[5:14,1], uiw = abs(coef(model)[5:14]-confint(model)[5:14,2]), add = T)


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
require(pROC)
set.seed(10)

subset
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