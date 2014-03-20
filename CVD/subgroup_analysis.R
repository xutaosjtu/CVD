## subgroup analysis
require(survival)

## S4
association_analysis = function(subset = which(S4$prev_mi==0)){
  rst = NULL;
  for (m in candidates){
    S4$metabolite = scale(log(S4[, m]))
    model = coxph(Surv(mi_time, inz_mi) ~ metabolite 
                   + scale(ltalteru) #+ as.factor(lcsex)
                   + scale(ltbmi)## model 1
                   + as.factor(my.diab)  ##model 2
                   + as.factor(my.alkkon) + as.factor(my.cigreg)
                   + scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
                   + scale(log(lh_crp))##model 4
#                   + as.factor(ltmstati)
#                   + as.factor(ltantihy)
#                   + as.factor(my.physical)
                  ,subset = subset,#
                  data = S4)
    rst = rbind(rst, summary(model)$coefficients[1,])
  }
  table(model$y[,2])
  rst = data.frame(rst,FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
  rst$lower = exp(rst$coef - 1.96*rst$se.coef.)
  rst$upper = exp(rst$coef+1.96*rst$se.coef.)
  rownames(rst) = candidates
  return(rst)
}

require(metafor)

heterogeneity = function(x,y){
  hetero = NULL;
  for(i in 1:length(candidates)){
    estimates = rbind(x[i,], y[i,])
    testrst = rma(yi = estimates$coef, sei= estimates$se.coef.,
                  method = "FE", 
                  measure="RR")
    summary = c(testrst$QE, testrst$QEp)
    hetero = rbind(summary, hetero)
  }
  colnames(hetero) = c("Heterogenity", "Hetero Pvalue")
  return(hetero)
}


## subgroup of age>65
asso.over60 = association_analysis(subset = S4$ltalter>=60&S4$prev_mi==0)
asso.below60 = association_analysis(subset = S4$ltalter<60&S4$prev_mi==0)
hetero = NULL
for(i in 1:length(candidates)){
  x = rbind(asso.over60[i,], asso.below60[i,])
  testrst = rma(yi = x$coef, sei= x$se.coef.,
                method = "FE", 
                measure="RR")
  summary = c(testrst$QE, testrst$QEp)
  hetero = rbind(summary, hetero)
}
colnames(hetero) = c("Heterogenity", "Hetero Pvalue")
write.csv(cbind(asso.over60, asso.below60, hetero), file = "metabolites_MI survival analysis_S4_Age subgroup.csv")

## subgroup of gender
asso.men = association_analysis(subset = S4$lcsex==1 & S4$prev_mi==0)
asso.women = association_analysis(subset = S4$lcsex==2 & S4$prev_mi==0)
hetero = heterogeneity(x=asso.men, y = asso.women)
write.csv(cbind(asso.men, asso.women, hetero), file = "metabolites_MI survival analysis_S4_sex subgroup.csv")

## subgroup of diseases
asso.diab = association_analysis(subset = S4$my.diab==1&S4$prev_mi==0) ## diabetes
asso.hyperten = association_analysis(subset = S4$lthyact==1&S4$prev_mi==0) ## hypertension
asso.none = association_analysis(subset = S4$lthyact==2&S4$my.diab==0&S4$prev_mi==0) 
write.csv(cbind(asso.diab, asso.hyperten, asso.none), file = "metabolites_MI survival analysis_S4_hypertension subgroup.csv")

## subgroup of smokers
asso.CS = association_analysis(subset = S4$my.cigreg==2 & S4$prev_mi==0) ## current smokers
asso.FS = association_analysis(subset = S4$my.cigreg==1 & S4$prev_mi==0) ## former smokers
asso.NS = association_analysis(subset = S4$my.cigreg==0 & S4$prev_mi==0) ## never smokers

write.csv(cbind(asso.CS, asso.FS, asso.NS), file = "metabolites_MI survival analysis_S4_smoker subgroup.csv")


## S2
S2.sub = S2[which(S2$subcoho==1|S2$inz_mi02==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Kynurenine_Trp","zz_nr","mi_time02", "inz_mi02","subcoho","prev_mi02","ctantihy", "ctmstati", "cl_ldla", "ctcigreg", "my.physical")]
#S2.sub = na.omit(S2.sub)
S2.sub = S2.sub[which(S2.sub$prev_mi02!=1),]
S2.sub$mi_time.start = 0
S2.sub$mi_time.start[which(S2.sub$inz_mi02==1)] = S2.sub$mi_time02[which(S2.sub$inz_mi02==1)]-1
S2.sub$mi_time.end = S2.sub$mi_time02

#cases in subcohort as control
#S2.sub$atcontrol = S2.sub$inz_mi
subcohort.cases = S2.sub[which(S2.sub$subcoho==1&S2.sub$inz_mi02==1),]
#subcohort.cases$atcontrol = 0
subcohort.cases$inz_mi02=0
subcohort.cases$mi_time.start = 0
subcohort.cases$mi_time.end = subcohort.cases$mi_time.end-1

S2.sub = rbind(S2.sub, subcohort.cases)

weight=rep(1,nrow(S2.sub))
weight[which(S2.sub$subcoho==1&S2.sub$inz_mi==0&S2.sub$ccsex==1)]= 1545/306
weight[which(S2.sub$subcoho==1&S2.sub$inz_mi==0&S2.sub$ccsex==2)]= 1699/289
#weight[which(S2.sub$subcoho==1 & S2.sub$inz_mi==0)]= 3244/813
#weight[which(S2.sub$subcoho==0 & S2.sub$inz_mi==1)]= (384-92)/87
#weight[which(S2.sub$subcoho==1 & S2.sub$inz_mi==1)] = 47/31
S2.sub$weight = weight

association_analysis = function(subset = 1:nrow(S2.sub)){
  rst = NULL
  for (m in candidates){
    S2.sub$metabolite = scale(log(S2.sub[, m]))
    model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ metabolite  
                  + scale(ctalteru) #+ as.factor(ccsex)## model 1
                  + scale(ctbmi) + as.factor(my.diab)  ##model 2
                  + scale(ctsysmm) #+ as.factor(my.cigreg)
                  + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
                  + scale(log(cl_crp))  ##model 4
                  #+ as.factor(ctmstati)
                  #+ cluster(as.factor(zz_nr))
                  #+ as.factor(my.physical)
                  ,data = S2.sub
                  ,subset = subset
                  ,weights = weight
                  ,method = "breslow"
    )
    rst = rbind(rst, summary(model)$coefficients[1,])
  }
  rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
  rst$lower = exp(rst$coef - 1.96*rst$se.coef.)
  rst$upper = exp(rst$coef+1.96*rst$se.coef.)
  rownames(rst) = candidates
  return(rst)
}

## Age <65
asso.over60 = association_analysis(subset = S2.sub$ctalteru>=60)
asso.below60 = association_analysis(subset = S2.sub$ctalteru<60)
hetero = heterogeneity(x = asso.over60, y = asso.below60)
write.csv(cbind(asso.over60, asso.below60, hetero), file = "metabolites_MI survival analysis_S2_Age subgroup.csv")

## subgroup of gender
asso.men = association_analysis(subset = S2.sub$ccsex==1)
asso.women = association_analysis(subset = S2.sub$ccsex==2)
hetero = heterogeneity(x = asso.men, y = asso.women)
write.csv(cbind(asso.men, asso.women, hetero), file = "metabolites_MI survival analysis_S2_sex subgroup.csv")

## subgroup of diseases
asso.diab = association_analysis(subset = S2.sub$my.diab==1) ## diabetes
asso.hyperten = association_analysis(subset = S2.sub$cthyact==1) ## hypertension
asso.none = association_analysis(subset = S2.sub$lthyact==2&S2.sub$my.diab==0) 
write.csv(cbind(asso.diab, asso.hyperten, asso.none), file = "metabolites_MI survival analysis_S4_hypertension subgroup.csv")

## subgroup of smokers
asso.CS = association_analysis(subset = S2.sub$my.cigreg==2) ## current smokers
asso.FS = association_analysis(subset = S2.sub$my.cigreg==1) ## former smokers
asso.NS = association_analysis(subset = S2.sub$my.cigreg==0) ## never smokers
write.csv(cbind(asso.CS, asso.FS, asso.NS), file = "metabolites_MI survival analysis_S2_smoker subgroup.csv")
