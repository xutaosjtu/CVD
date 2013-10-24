require(survival)
require(PredictABEL)
require(pROC)
require(boot)

S2$Arg.Trp = S2$Arg/S2$Trp

data.S2 = subset(S2, subcoho==1 | inz_mi02==1)
data.S2 = subset(data.S2, prev_mi02==0)
data.S2 = data.S2[which(data.S2$prev_mi02!=1),]

#
data.S2$prev_mi = data.S2$prev_mi02
data.S2$inz_mi = data.S2$inz_mi02
data.S2$mi_time = data.S2$mi_time02
#

data.S2$mi_time.start = 0
data.S2$mi_time.start[which(data.S2$inz_mi==1)] = data.S2$mi_time[which(data.S2$inz_mi==1)]-1
data.S2$mi_time.end = data.S2$mi_time

#cases in subcohort as control
subcohort.cases = data.S2[which(data.S2$subcoho==1&data.S2$inz_mi==1),]
#subcohort.cases$atcontrol = 0
subcohort.cases$inz_mi=0
subcohort.cases$mi_time.start = 0
subcohort.cases$mi_time.end = subcohort.cases$mi_time.end-1

data.S2 = rbind(data.S2, subcohort.cases)

weight=rep(1,nrow(data.S2))
weight[which(data.S2$subcoho==1&data.S2$inz_mi==0&data.S2$ccsex==1)]= 1545/306
weight[which(data.S2$subcoho==1&data.S2$inz_mi==0&data.S2$ccsex==2)]= 1699/289
#weight[which(tmp$subcoho==0 & tmp$inz_mi==1)]= (384-92)/87
#weight[which(tmp$subcoho==1 & tmp$inz_mi==1)] = 92/59
data.S2$weight = weight

#data.S2 = na.omit(data.S2[,c(clinical.S2, S2_valid_measures)])

clinical.S2 = c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi")

#################  Framingham score	##########################
framingham<-function(x, method="linear"){
  #print(x["lcsex"])
  if(x["ccsex"] == 1){ beta = c(3.06117, 1.12370, -0.93263, 1.99881, 1.93303, 0.65451, 0.57367, 23.9802); S0 = 0.88936}
  else {beta = c(2.32888, 1.20904, -0.70833, 2.76157, 2.82263, 0.52873, 0.69154, 26.1931); S0 = 0.95012}
  X = rep(0, 6)
  X[1] = beta[1] * log(x$ctalteru)#age
  X[2] = beta[2] * log(x$cl_chola)#total cholesterol
  X[3] = beta[3] * log(x$cl_hdla)#HDL
  if(is.na(x$ctantihy)) { X[4]=NA }
  else if(!(x$ctantihy ==1)) X[4]= beta[4]*log(x$ctsysmm) #systolic blood pressure if un-treated
  else  X[4] = beta[5]*log(x$ctsysmm) #systolic blood pressure if treated
  X[5] = beta[6]* as.numeric(x$ctcigreg==1|x$ctcigreg==2)#smoking
  X[6] = beta[7] * (as.numeric(x$my.diab))
  if(method=="score") return(1-S0^exp(sum(X) - beta[8]))
  else if(method == "linear") return(sum(X) - beta[8])
}

for(i in 1:nrow(data.S2)){
  data.S2$framingham.linear[i] = framingham(data.S2[i,], method="linear") 
}

###	men and women separated, using the framingham score
data = data.frame(
  start = data.S2$mi_time.start, end = data.S2$mi_time.end,
  event = data.S2$inz_mi,  # time and events
  framingham.linear = data.S2$framingham.linear,
  scale(log(as.matrix(data.S2[, S2_valid_measures]))),
  weight = data.S2$weight
)
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

# loglik = 0
model = coxph(Surv(start, end, event) ~ .,
              weights = data$weight,
                #subset = which(data.S2$ccsex ==i),
                data[ ,c( "start", "end", "event", metabo.selected,"framingham.linear")])
prediction = predict(model, type="risk")
fits[[2]] = roc(model$y[,3],prediction)
#  print(data.frame("lcsex"= i, "basline" = sort(survfit(model)$surv)[1]))
#  prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
#  loglik = model$loglik[2] + loglik




# 
# ##likelihood ratio test
# logliks = list()
# logliks[[1]] = loglik
# logliks[[2]] = loglik
# pchisq(-2*(logliks[[1]] - logliks[[2]]), df=1) 

theta.fit <- function(x, y, weight,...) 
{
  #d = data.frame(y, x)
  #print(dim(d))
  #print(colnames(d))
  coxph(y ~  ., weights = weight, data = as.data.frame(x))
}
theta.predict <- function(fit, x)
{
  #if(is.null(dim(x))) x=t(x)
  #dim(x)
  #print(colnames(x))
  value=predict(fit,newdata= as.data.frame(x) , type="risk")
  S0=min(survfit(fit)$surv)
  value = 1- S0^value
  return(value)
}
crossval.cox = function (x, y, theta.fit, theta.predict, weight, ..., ngroup = n) 
{
  call <- match.call()
  #x <- as.matrix(x)
  if(mode(x)=="numeric")  n = length(x)
  else n = nrow(x)
  ngroup <- trunc(ngroup)
  if (ngroup < 2) {
    stop("ngroup should be greater than or equal to 2")
  }
  if (ngroup > n) {
    stop("ngroup should be less than or equal to the number of observations")
  }
  if (ngroup == n) {
    groups <- 1:n
    leave.out <- 1
  }
  if (ngroup < n) {
    leave.out <- trunc(n/ngroup)
    o <- sample(1:n)
    groups <- vector("list", ngroup)
    for (j in 1:(ngroup - 1)) {
      jj <- (1 + (j - 1) * leave.out)
      groups[[j]] <- (o[jj:(jj + leave.out - 1)])
    }
    groups[[ngroup]] <- o[(1 + (ngroup - 1) * leave.out):n]
  }
  u <- NULL
  cv.fit <- rep(NA, n)
  for (j in 1:ngroup) {
    s = rep(T, n); s[groups[[j]]]=F
    u <- theta.fit(subset(x, subset = s), y[-groups[[j]],], weight[-groups[[j]]])
    cv.fit[groups[[j]]] <- theta.predict(u, subset(x, subset = !s))
  }
  if (leave.out == 1) 
    groups <- NULL
  return(list(cv.fit = cv.fit, ngroup = ngroup, leave.out = leave.out, 
              groups = groups, call = call))
}


prediction = rep(NA, nrow(data))
names(prediction) = rownames(data)
## combine men and women
subset = setdiff(which(!is.na(data.S2$ccsex)), na.index)# & S4$ltmstati==2, "framingham.linear"
pred.cv = crossval.cox(x = data[subset, c(metabo.selected, "framingham.linear")], y= Surv(data$start[subset], data$end[subset], data$event[subset]), theta.fit,theta.predict, weight=data$weight[subset], ngroup = length(subset))
prediction[subset] = pred.cv$cv.fit 

metabo.selected,

##men
subset = setdiff(which(data.S2$ccsex ==1), na.index)# & S4$ltmstati==2
pred.cv = crossval.cox(x = data[subset, c("framingham.linear")], y= Surv(data$start[subset], data$end[subset], data$event[subset]), theta.fit, theta.predict, weight=data$weight[subset], ngroup = length(subset))
prediction[subset] = pred.cv$cv.fit 
##women
subset = setdiff(which(data.S2$ccsex ==2), na.index)# & S4$ltmstati==2, "framingham.linear"
pred.cv = crossval.cox(x = data[subset, c("framingham.linear")], y= Surv(data$start[subset], data$end[subset], data$event[subset]), theta.fit, theta.predict, weight=data$weight[subset], ngroup = length(subset))
prediction[subset] = pred.cv$cv.fit 

fits = list()
fits[[1]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
fits[[2]] = roc (data$event[which(!is.na(prediction))], data.S2$framingham.score[which(!is.na(prediction))], ci = T)
#fits[[2]] = roc (data$event[which(!is.na(prediction))], data.S2$framingham.score[which(!is.na(prediction))], ci = T)


#calculate delta AUC
fits.test =  roc.test(fits[[1]], fits[[2]])
deltaAUC <- function(fits.test){
  delta = abs(fits.test$estimate[1] - fits.test$estimate[2]) 
  se = delta / fits.test$statistic
  return(data.frame("deltaAUC" = delta, "lower" = delta - 1.96*se, "upper" = delta + 1.96*se))
}
deltaAUC(fits.test)

#calculate NRI and IDI
reclassification(data[which(!is.na(prediction)), ], cOutcome = 3, fits[[2]]$predictor, fits[[1]]$predictor, cutoff = c(0, 0.10, 0.2, 1))

##plot
plot(fits[[2]], main = "Framingham Score", lty=2,legacy.axes=T)
plot(fits[[1]], col="red",add = T,lty=1,legacy.axes=T)
auc.ref = paste(round(fits[[2]]$auc[1],2),"(",round(fits[[2]]$ci[1],2),",", round(fits[[2]]$ci[3],2),")", sep="")
auc.full = paste(round(fits[[2]]$auc[1],2),"(",round(fits[[1]]$ci[1],2),",", round(fits[[1]]$ci[3],2),")", sep="")
legend(0.8, 0.2, 
       legend = paste(c("Ref","Plus metabo"), c(auc.ref, auc.full)), 
       col = c("black","red"), lty = c(2:1)
)  

###############################################################
######  combined men and women, using reference model 4 #######
###############################################################
data = data.frame(
  start = data.S2$mi_time.start, end = data.S2$mi_time.end, event = data.S2$inz_mi,  # time and events
  scale(data.S2$ctalteru),  as.factor(data.S2$ccsex), ##model 1
  scale(data.S2$ctbmi),as.factor(data.S2$my.diab), ##model 2
  scale(data.S2$ctsysmm),  as.factor(data.S2$my.cigreg), as.factor(data.S2$my.alkkon), scale(data.S2$cl_hdla), scale(data.S2$cl_chola), ##model 3
  scale(data.S2$cl_crp), ##model 4,
  scale(log(as.matrix(data.S2[, S2_valid_measures]))),
  weight=data.S2$weight
)
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

clinical = c("ctalteru",  "ccsex", "ctbmi", "my.diab", "ctsysmm", "my.cigreg", "my.alkkon", "cl_hdla", "cl_chola", "cl_crp")
colnames(data)[4:13] = clinical
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

ref = list()
ref[[1]] = clinical[1:2]
ref[[2]] = clinical[1:4]
ref[[3]] = clinical[1:9]
ref[[4]] = clinical
# #estimation without cross-validation
# prediction = rep(NA, dim(data)[1])
# names(prediction) = rownames(data)
# loglik = 0
# 
# model = coxph(Surv(data$time, data$event) ~ .
#               ,data[, c( ref[[4]])])
# prediction[dimnames(model$y)[[1]]] =  1 - sort(survfit(model)$surv)[1] ^predict(model, type = "risk")
# loglik = model$loglik[2]
# sort(survfit(model)$surv)[1]
# metabo.selected3,, "ltmstati"
# "Arg"
# 
# logliks = list()
# logliks[[1]] = loglik
# logliks[[2]] = loglik
# pchisq(-2*(logliks[[1]] - logliks[[2]]), df=4)



for(i in 1:4){
  prediction = rep(NA, dim(data)[1])
  names(prediction) = rownames(data)
  subset = setdiff(which(data.S2$prev_mi == 0), na.index)#& S4$ltmstati!=1, "ltmstati"
  pred = crossval.cox(x = data[subset, c(ref[[i]])], y= Surv(data$start[subset], data$end[subset], data$event[subset]), theta.fit, theta.predict, weight=data$weight[subset], ngroup = length(subset))
  prediction[subset] = pred$cv.fit
  fits[[i]] = roc (data$event[which(!is.na(prediction))], prediction[which(!is.na(prediction))], ci = T)
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

#calculate NRI and IDI
require(PredictABEL)
for(i in 1:4){
  print(paste("Evaluation of model", i))
  reclassification(data[which(!is.na(prediction)), ], cOutcome =3, fits.ref[[i]]$predictor, fits.full[[i]]$predictor, cutoff = c(0, 0.1, 0.2, 1))
}

##plots
pdf("added predictive value in S2 case cohort.pdf",width=8, height=8)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(fits.ref[[i]], main = paste("Model",i), lty=2,legacy.axes=T)
  plot(fits.full[[i]], col="red",add = T,lty=1,legacy.axes=T)
  auc.ref = paste(round(fits.ref[[i]]$auc[1],2),"(",round(fits.ref[[i]]$ci[1],2),",", round(fits.ref[[i]]$ci[3],2),")", sep="")
  auc.full = paste(round(fits.full[[i]]$auc[1],2),"(",round(fits.full[[i]]$ci[1],2),",", round(fits.full[[i]]$ci[3],2),")", sep="")
   legend(0.8, 0.2, 
          legend = paste(c("Ref","Plus metabo"), c(auc.ref, auc.full)), 
          col = c("black","red"), lty = c(2:1)
   )  
}
dev.off()









##estimates of the confounders and prediction model evaluation by sample spliting
## sample spliting
training = sample(1:692, 461)
test = setdiff(c(1:692), training)
## model fitting
model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ 
                scale(log(Arg)) + 
                scale(log(lysoPC_a_C17_0)) + 
                scale(log(Trp)) + 
                scale(log(PC_aa_C32_2))+
                #                 scale(log(PC_aa_C36_3))+
                #                 scale(log(lysoPC_a_C18_2)) + 
                #                 scale(log(SM_C24_1))+
                scale(ctalteru) + as.factor(ccsex)
              + scale(ctbmi) + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
              + scale((cl_crp))  ##model 
              #+ cluster(as.factor(zz_nr))
              ,data = S2.sub[training,]#[which(S2.sub$ctmstati!=1)]
              ,weights = weight
              ,method = "breslow"
)
#write.csv(cbind(summary(model)$coef, confint(model)), file = "estimates of confounders plus original four metabolites in S2_model 4_remove statine.csv")
model.base = update(model, .~. - scale(log(Arg))-scale(log(Trp))-scale(log(lysoPC_a_C17_0))-scale(log(PC_aa_C32_2)))
pred = predict(model.base, S2.sub[test,], type="risk")
S0=min(survfit(model.base)$surv)
pred = 1- S0^pred
fits[[1]] = roc(S2.sub[test, 'inz_mi02'], pred, ci= T)
pred = predict(model, S2.sub[test,], type="risk")
S0=min(survfit(model)$surv)
pred = 1- S0^pred
fits[[2]] = roc(S2.sub[test, 'inz_mi02'], pred, ci= T)

reclassification(S2.sub[test[which(!is.na(pred))], ], cOutcome = 146, fits[[1]]$predictor, fits[[2]]$predictor, cutoff = c(0, 0.10, 0.2, 0.5,1))

for(i in 1:nrow(S2.sub)){
  S2.sub$framingham.score[i] = framingham(S2.sub[i,], method="score") 
}

model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ 
                scale(log(Arg)) + 
                scale(log(lysoPC_a_C17_0)) + 
                scale(log(Trp)) + 
                scale(log(PC_aa_C32_2))+
                framingham.linear
              ,data = S2.sub[training,]#[which(S2.sub$ctmstati!=1)]
              ,weights = weight
              ,method = "breslow"
)
model.base = update(model, .~. - scale(log(Arg))-scale(log(Trp))-scale(log(lysoPC_a_C17_0))-scale(log(PC_aa_C32_2)))
pred = predict(model.base, S2.sub[test,], type="risk")
S0=min(survfit(model.base)$surv)
pred = 1- S0^pred
fits[[1]] = roc(S2.sub[test, 'inz_mi02'], pred, ci = T)
pred = predict(model, S2.sub[test,], type="risk")
S0=min(survfit(model)$surv)
pred = 1- S0^pred
fits[[2]] = roc(S2.sub[test, 'inz_mi02'], pred, ci = T)

reclassification(S2.sub[test[which(!is.na(pred))], ], cOutcome = 146, fits[[1]]$predictor, fits[[2]]$predictor, cutoff = c(0, 0.10, 0.2, 0.5,1))
