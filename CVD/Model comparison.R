# TODO: Add comment
# 
# Author: tao.xu
###############################################################################
require(survival)
require(PredictABEL)
require(pROC)
require(boot)

S4$my.sysmm.untreat = log(S4$ltsysmm)
S4$my.sysmm.untreat[which((S4$ltmbbl ==1 | S4$ltmace ==1 | S4$ltmata ==1))] = 0
S4$my.sysmm.treat = log(S4$ltsysmm)
S4$my.sysmm.treat[which(!(S4$ltmbbl ==1 | S4$ltmace ==1 | S4$ltmata ==1))] = 0
data = subset(S4, prev_mi==0& !is.na(inz_mi))
data[, metabo.selected3] =  log(S4[, metabo.selected3])

#################	Framingham score	##########################
framingham<-function(x, method="linear"){
  #print(x["lcsex"])
  if(x["lcsex"] == 1){ beta = c(3.06117, 1.12370, -0.93263, 1.99881, 1.93303, 0.65451, 0.57367, 23.9802); S0 = 0.88936}
  else {beta = c(2.32888, 1.20904, -0.70833, 2.76157, 2.82263, 0.52873, 0.69154, 26.1931); S0 = 0.95012}
  X = rep(0, 6)
  X[1] = beta[1] * log(x$ltalteru)#age
  X[2] = beta[2] * log(x$ll_chola)#total cholesterol
  X[3] = beta[3] * log(x$ll_hdla)#HDL
  if(is.na(x$ltantihy)) { X[4]=NA }
  else if(!(x$ltantihy ==1)) X[4]= beta[4]*log(x$ltsysmm) #systolic blood pressure if un-treated
  else  X[4] = beta[5]*log(x$ltsysmm) #systolic blood pressure if treated
  X[5] = beta[6]* as.numeric(x$ltcigreg==1|x$ltcigreg==2)#smoking
  X[6] = beta[7] * (as.numeric(x$my.diab))
  if(method=="score") return(1-S0^exp(sum(X) - beta[8]))
  else if(method == "linear") return(sum(X) - beta[8])
}

for(i in 1:nrow(S4)){
	S4$framingham[i] = framingham(S4[i,], method = "score") 
}

for(i in 1:nrow(S4)){
  S4$framingham.linear[i] = framingham(S4[i,], method = "linear") 
}

###	men and women separated, using the framingham score
data = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		framingham.linear = S4$framingham.linear,
		scale(log(as.matrix(S4[, S4_valid_measures]))) , 
		lh_crp = scale(S4$lh_crp)
)
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))
prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)

loglik = 0
for (i in 1:2){
	model = coxph(Surv(time, event) ~ .,
			subset = which(S4$prev_mi == 0&S4$lcsex ==i),#$ltmstati==2
			data[ ,c( "time", "event", metabo.selected3, "framingham.linear")])
	print(data.frame("lcsex"= i, "basline" = sort(survfit(model)$surv)[1]))
	prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
	loglik = model$loglik[2] + loglik
}

##likelihood ratio test
logliks = list()
logliks[[1]] = loglik
logliks[[2]] = loglik
pchisq(-2*(logliks[[1]] - logliks[[2]]), df=1) 

##men
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$lcsex ==1), na.index)# & S4$ltmstati==2
pred.cv = crossval.cox(x = data[subset, c(metabo.selected3, "lh_crp", "framingham.linear")], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] = 1-0.6742013^ pred.cv$cv.fit 
##women
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$lcsex ==2), na.index)# & S4$ltmstati==2
pred.cv = crossval.cox(x = data[subset, c(metabo.selected3,"lh_crp", "framingham.linear")], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] =1- 0.9846102 ^ pred.cv$cv.fit 

fits = list()
fits[[1]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
#fits[[2]] = roc (data$inz_mi[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
fits[[2]] = roc (data$event[which(!is.na(prediction))], S4$framingham[which(!is.na(prediction))], ci = T)

## ROC
pdf("ROC_Framingham.pdf")
plot(fits.framingham[[1]]) 
plot(fits.framingham[[2]], col = "green", add = T)
plot(fits.framingham_crp[[1]], col = "blue", add = T)
plot(fits.framingham_crp[[3]], col = "red", add = T)
dev.off()

#calculate delta AUC
fits.test =  roc.test(fits[[1]], fits[[2]])
deltaAUC <- function(fits.test){
	delta = abs(fits.test$estimate[1] - fits.test$estimate[2]) 
	se = delta / fits.test$statistic
	return(data.frame("deltaAUC" = delta, "lower" = delta - 1.96*se, "upper" = delta + 1.96*se))
}
deltaAUC(fits.test)

#calculate NRI and IDI
fits.recl = reclassification(data[which(!is.na(prediction)), ], cOutcome = 2, fits[[2]]$predictor, fits[[1]]$predictor, cutoff = c(0, 0.1, 0.2,  1))


###	men and women separated, using the framingham model
data = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		log(S4$ltalteru), #age
		as.factor(S4$lcsex), #sex
		log(S4$ll_chola),	#total cholesterol
		log(S4$ll_hdla),	#HDL cholesterol
		S4$my.sysmm.untreat,#systolic blood pressure untreated
		S4$my.sysmm.treat,	#systolic blood pressure treated    
		as.numeric(S4$ltcigreg==1|S4$ltcigreg==2),	#smoking 
		(as.numeric(S4$my.diab)-1),	#diabetes
		scale(log(S4[, metabo.asso]))#, scale(S4[, metabo.ratio.asso])
)
clinical = c("age", "sex",  "ll_chola", "ll_hdla", "my.sysmm.untreat", "my.sysmm.treat", "smoking", "diabetes")
colnames(data)[3:10] = c("age", "sex",  "ll_chola", "ll_hdla", "my.sysmm.untreat", "my.sysmm.treat", "smoking", "diabetes")
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))
prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)

##without cross-validation
loglik = 0
for( i in 1:2){
	model = coxph(Surv(data$time, data$event) ~ .,
		subset = which(S4$prev_mi == 0 & S4$lcsex ==i & S4$ltmstati==2),
		data[, c(metabo.selected3, clinical[-2])])
	print(data.frame("lcsex"= i, "basline" = sort(survfit(model)$surv)[1]))
	prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
	loglik = model$loglik[2] + loglik
}
metabo.selected3,
##likelihood ratio test
logliks = list()
logliks[[1]] = loglik
logliks[[2]] = loglik
pchisq(-2*(logliks[[1]] - logliks[[2]]), df=1) 

#with cross-validation
#men
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$lcsex ==1 ), na.index)#& S4$ltmstati==2
pred.cv = crossval.cox(x = data[subset, c(clinical[-2], metabo.selected3)], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] = 1 -  0.6187485 ^ pred.cv$cv.fit #metabolite+reference:0.6187485; reference: 0.6882019
#exclude statin ref: 0.6819058, ref+metabo: 0.6390437
#women
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$lcsex ==2 ), na.index)#& S4$ltmstati==2
pred.cv = crossval.cox(x = data[subset, c(clinical[-2], metabo.selected3)], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] =1-   0.9888222 ^ pred.cv$cv.fit #metabolite+reference:0.9888222; reference: 0.9735755
#exclude statin ref: 0.9729412, ref + metabo: 0.9869769
, metabo.selected3


fits = list()
fits[[1]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
fits[[2]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)

#calculate delta AUC
fits.test =  roc.test(fits[[1]], fits[[2]])
deltaAUC <- function(fits.test){
	delta = abs(fits.test$estimate[1] - fits.test$estimate[2]) 
	se = delta / fits.test$statistic
	return(data.frame("deltaAUC" = delta, "lower" = delta - 1.96*se, "upper" = delta + 1.96*se))
}
deltaAUC(fits.test)

#calculate NRI and IDI
reclassification(data[which(!is.na(prediction)), ], cOutcome = 2, fits[[2]]$predictor, fits[[1]]$predictor, cutoff = c(0, 0.1, 0.2, 1))

##########################	reference model used to selecte metabolites	########################
###	men and women separated, using reference model 4
prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)
loglik = 0
#without cross-validation
for( i in 1:2){
	model = coxph(Surv(mi_time, inz_mi) ~ scale(Arg.Trp) + scale(log(PC_aa_C32_2)) + scale(log(lysoPC_a_C17_0))+ #scale(log(Trp))+scale(log(Arg)) + scale(Arg.Trp) + scale(log(Trp))+scale(log(Arg)) 
					log(ltalteru) + as.factor(lcsex) +log(ltbmi)## model 1
					+ my.diab ##model 2
					+log(ltsysmm) + my.cigreg  + my.alkkon + log(ll_chola) + log(ll_hdla) ##model 3
					+ scale(lh_crp)  ##model 4 + my.physical
			,subset = which(S4$prev_mi == 0 & S4$ltmstati!=1),#&S4$lcsex ==i, S4$ltmstati!=1
			S4)
	print(data.frame("lcsex"= i, "basline" = sort(survfit(model)$surv)[1]))
	prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
	loglik = model$loglik[2] + loglik
}

coxph(Surv(time,event)~., data = tmp2[, c("time", "event", names(coef(model.penal.opt$fullfit)))])

rst = vector(mode = "numeric", length = dim(S4)[1])
name = ""
for(i in 1:(length(S4_valid_measures)-1)){
	tmp = log(S4[, S4_valid_measures[i]]) - log(S4[, S4_valid_measures[(i+1):length(S4_valid_measures)]])
	name = c(name, paste(S4_valid_measures[i], S4_valid_measures[(i+1):length(S4_valid_measures)], sep = "."))
	rst = cbind(rst, tmp)
}


##likelihood ratio test
logliks = list()
logliks[[1]] = loglik
logliks[[2]] = loglik
1-pchisq(2*(logliks[[1]] - logliks[[2]]), df=4) 

###with cross validation
data = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		S4$ltalteru, as.factor(S4$lcsex), S4$ltbmi, ##model 1
		S4$my.diab, ##model 2
		S4$ltsysmm, S4$ll_hdla, S4$ll_chola, as.factor(S4$my.cigreg), S4$my.alkkon,##model 3
		log(S4$lh_crp), ##model 4
		log(as.matrix(S4[, metabo.asso]))
)
clinical = c("age",  "sex", "ltbmi","diabetes", "ltsysmm", "ll_hdla", "ll_chola", "smoking", "alkkon", "lh_crp")
colnames(data)[3:12] = clinical#, "total2HDL"
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

#men
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) ), na.index) #& S4$lcsex ==1
pred = crossval.cox(x = data[subset, c(metabo.selected3,clinical[-3])], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] = 1- 0.5676883 ^ pred$cv.fit
#women
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$lcsex ==2), na.index)
pred = crossval.cox(x = data[subset, c(metabo.selected3, clinical[-3])], y= Surv(tmp$time[subset], tmp$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
prediction[subset] = 1- 0.989219 ^ pred$cv.fit

fits = list()
fits[[1]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
fits[[2]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)

#calculate delta AUC
fits.test =  roc.test(fits[[1]], fits[[2]])
deltaAUC <- function(fits.test){
	delta = abs(fits.test$estimate[1] - fits.test$estimate[2]) 
	se = delta / fits.test$statistic
	return(data.frame("deltaAUC" = delta, "lower" = delta - 1.96*se, "upper" = delta + 1.96*se))
}
deltaAUC(fits.test)

#calculate NRI and IDI
#NRI(category): sum(v[i])/count(EVENTS)-sum(v[j])/count(NONEVENTS)
#NRI(continous): sum(pnew[i]-pold[i]) - sum(pnew[j] - pold[j])
#IDI = (IS[new] - IS [old]) - (IP[new] - IP[old])
reclassification(data[which(!is.na(prediction)), ], cOutcome = 2, fits[[2]]$predictor, fits[[1]]$predictor, cutoff = c(0, 0.03, 0.08, 0.15, 1))

#######################################################################
###	combined men and women, using reference models
prediction = rep(NA, dim(data)[1])
names(prediction) = rownames(data)
loglik = 0
data = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		scale(S4$ltalteru),  as.factor(S4$lcsex), scale(S4$ltbmi),##model 1
		S4$my.diab, ##model 2
		scale(S4$ltsysmm),  as.factor(S4$my.cigreg), as.factor(S4$my.alkkon), scale(S4$ll_hdla), scale(S4$ll_chola), ##model 3
		scale(log(S4$lh_crp)), ##model 4,
		'ltmstati'=as.factor(S4$ltmstati), ##adding statin as a covariate
		scale(log(S4[, S4_valid_measures]))#, scale(S4[,metabo.ratio.asso])
)
clinical = c("age", "sex", "ltbmi", "diabetes", "ltsysmm", "smoking", "alkkon", "ll_hdla", "ll_chola", "lh_crp")
colnames(data)[3:12] = clinical#, "total2HDL"
#na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

ref = list()
ref[[1]] = clinical[1:2]
ref[[2]] = clinical[1:4]
ref[[3]] = clinical[1:9]
ref[[4]] = clinical
#estimation without cross-validation
model = coxph(Surv(data$time, data$event) ~ .
		,subset = which(S4$prev_mi == 0),#&S4$ltmstati !=1
		data[, c( ref[[4]])])
prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
loglik = model$loglik[2]
sort(survfit(model)$surv)[1]
metabo.selected3,, "ltmstati"
"Arg"

logliks = list()
logliks[[1]] = loglik
logliks[[2]] = loglik
pchisq(-2*(logliks[[1]] - logliks[[2]]), df=4)


fits=list()
for(i in 1:4){
	prediction = rep(NA, dim(data)[1])
	names(prediction) = rownames(data)
	subset = which(S4$prev_mi == 0 & !is.na(S4$inz_mi))#, na.index)#& S4$ltmstati!=1, "ltmstati"
	pred = crossval.cox(x = data[subset, c(ref[[i]])], 
                      y= Surv(data$time[subset], data$event[subset]), 
                      theta.fit, theta.predict, ngroup = length(subset)
                      )
	prediction[subset] = 1-0.856785 ^ pred$cv.fit
	fits[[i]] = roc (data$event[which(!is.na(prediction))], 
                   prediction[which(!is.na(prediction))], ci = T)
}




deltaAUC <- function(fits.test){
	delta = abs(fits.test$estimate[1] - fits.test$estimate[2]) 
	se = delta / fits.test$statistic
	return(data.frame("deltaAUC" = delta, "lower" = delta - 1.96*se, "upper" = delta + 1.96*se))
}


#calculate delta AUC
for(i in 1:4){
	#print(fits[[i]]$ci[1:3])
	#fits[[5]]$ci[1:3]
	fits.test =  roc.test(fits.full[[i]], fits.ref[[i]])
	print(fits.test)
	print(deltaAUC(fits.test))
}


## plot ROC in the referenc and full model
pdf("added predictive value 3 markers.pdf")
plot(fits.full[[4]])
plot(fits.ref[[3]], lty = 2, add = T)
dev.off()

pdf("ROC curves_S4.pdf")
plot(fits.ref[[3]]) ## without metabolite and CRP
plot(fits.ref[[4]], add = T, lty = 2, col = "blue") ## with CRP, without metabolites
plot(fits.full[[4]], lty = 3, add = T, col = "red") ## with CRP and metabolites
plot(fits.full[[3]], lty = 4, add = T, col = "green") ## with metabolites without CRP
dev.off()

#calculate NRI and IDI
require(PredictABEL)
for(i in 1:4){
	print(paste("Evaluation of model", i))
	reclassification(data[which(!is.na(prediction)), ], cOutcome = 2, fits.ref[[i]]$predictor, fits.full[[i]]$predictor, cutoff = c(0, 0.1, 0.2, 1))
}

## adding sigle metabolite(ratio) into the model
for(i in 1:3){
	prediction = rep(NA, dim(data)[1])
	names(prediction) = rownames(data)
	subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi) & S4$ltmstati!=1), na.index)
	pred = crossval.cox(x = data[subset, c(metabo.selected3[i],ref[[4]])], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
	prediction[subset] = 1-0.8482743 ^ pred$cv.fit
	fits[[i]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
}

for(i in 1:3){
	fits.test =  roc.test(fits.ref4[[i]], fits.dif_ref[[4]])
	print(fits.test)
	print(deltaAUC(fits.test))
}

#calculate NRI and IDI
require(PredictABEL)
for(i in 1:3){
	print(paste("Evaluation of model", i))
	reclassification(data[which(!is.na(prediction)), ], cOutcome = 2, fits.dif_ref[[4]]$predictor, fits.ref4[[i]]$predictor, cutoff = c(0, 0.03, 0.08, 0.15, 1))
}

metabo.selected3[i],

<- function(metabolite, ref, data, subset){
	 model = coxph(Surv(data$time, data$event) ~ .
			 ,subset = which(S4$prev_mi == 0),
			 data[, c(metabolite, ref)])
	 baseSurv = sort(survfit(model)$surv)[1]
	 prediction = rep(NA, dim(data)[1])
	 names(prediction) = rownames(data)
	 #subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi)), na.index)
	 pred = crossval.cox(x = data[subset, c(metabo.selected3, clinical)], y= Surv(data$time[subset], data$event[subset]), theta.fit, theta.predict, ngroup = length(subset))
	 prediction[subset] = 1-   baseSurv ^ pred$cv.fit
}

fits = list()
fits[[1]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
fits[[2]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)


#subset =  which(S4$prev_mi == 0&S4$lcsex ==2)
#pred.3 = roc (data$inz_mi[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
#
#pred.framingham = roc (data$inz_mi[which(S4$prev_mi == 0)], data$framingham[which(S4$prev_mi == 0)], ci = T)
#
#base = data.frame(ltalteru=0, ll_chola =0, ll_hdla = 0, my.sysmm.untreat = 0, my.sysmm.treat = 0, ltcigreg = 0, my.diab = 0)

extractAIC(model)

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

#########plot of ROC
par(mfrow = c(2,2))
for (i in 1:4){
	plot(fits.dif_ref[[i]], lwd = 2, lty = 2, col = "blue", main = paste("Model", i))
	plot(fits.dif_full[[i]], lwd = 2, lty = 1, col = "red", add = T)
	ciref = ci(fits.dif_ref[[i]])
	cifull = ci(fits.dif_full[[i]])
	legend(0.6, 0.2, 
			legend = c(
					substitute(Reference (AUC: k (l,u)), 
							list(k = round(ciref[2],2), 
									l = round(ciref[1],2), 
									u = round(ciref[3],2))
					),
					substitute(Full (AUC: k (l,u)), 
							list(k = round(cifull[2],2),
									l = round(cifull[1],2),
									u = round(cifull[3],2))
					)
			), 
			col = c("blue","red"), lty = c(2,1))
}


###
model.female = coxph(Surv(mi_time, inz_mi) ~ #scale(log(Arg))
                     #+ as.factor(ltnuecht):scale(log(Trp))
                     + scale(log(lysoPC_a_C18_1)) 
                     + scale(log(PC_ae_C38_0)) + 
               scale(ltalteru) 
              + scale(ltbmi)## model 1
              + my.diab  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(ll_chola) + scale(ll_hdla) 
              + scale(lh_crp)  ##model 4
                    # + strata(ltnuecht)
              ,subset = which(S4$prev_mi==0 & S4$lcsex ==2 & S4$ltnuecht==1),
              data = S4)

model.male = coxph(Surv(mi_time, inz_mi) ~ #scale(log(Arg)) 
                   #+ as.factor(ltnuecht):scale(log(Trp))
                   + scale(log(lysoPC_a_C18_1)) 
                   + scale(log(PC_ae_C38_0)) + 
                      scale(ltalteru) 
                     + scale(ltbmi)## model 1
                     + my.diab  ##model 2
                     + scale(ltsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(ll_chola) + scale(ll_hdla)
                     + scale(lh_crp)  ##model 4
                   #+ strata(ltnuecht)
                     ,subset = which(S4$prev_mi==0 & S4$lcsex ==1 & S4$ltnuecht==1),
                     data = S4)
sort(survfit(model.female)$surv)[1]
sort(survfit(model.male)$surv)[1]

pred = S4$inz_mi
fit1 = roc(c(model.male$y[,2], model.female$y[,2]), c((1- 0.6832817^predict(model.male,type="risk")), (1-0.9560963^predict(model.female, type="risk"))), ci =T )

pred = S2$inz_mi;
male= which(S2$prev_mi==0 & S2$lcsex ==1);
female=which(S2$prev_mi==0 & S2$lcsex ==2)
pred[male] = 1- 0.6832817^predict(model.male,newdata=S2[male,],type="risk")
pred[female] = 1-0.9560963^predict(model.female,newdata=S2[female,], type="risk")
fit2 = roc(S2$inz_mi[c(male,female)], pred[c(male,female)], ci = T)

newdata=S2[female,],

scale(log(Arg))+ as.factor(ltnuecht):scale(log(Arg))
#+ as.factor(ltnuecht):scale(log(Trp))
+ scale(log(lysoPC_a_C18_1)) + as.factor(ltnuecht):scale(log(lysoPC_a_C18_1))
+ scale(log(PC_ae_C38_0)) + as.factor(ltnuecht):scale(log(PC_ae_C38_0)) + 

  scale(log(Arg))+ as.factor(ltnuecht):scale(log(Arg))
#+ as.factor(ltnuecht):scale(log(Trp))
+ scale(log(lysoPC_a_C18_1)) + as.factor(ltnuecht):scale(log(lysoPC_a_C18_1))
+ scale(log(PC_ae_C38_0)) + as.factor(ltnuecht):scale(log(PC_ae_C38_0)) + 
  