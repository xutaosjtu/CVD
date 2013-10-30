S2$Arg.Trp = S2$Arg/S2$Trp

data.S2 = subset(S2, subcoho==1 | inz_mi02==1)
data.S2 = subset(data.S2, prev_mi02==0)

clinical.S2 = c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi")

data.S2 = na.omit(data.S2[,c(clinical.S2, S2_valid_measures)])

# S2.sub = S2[which((S2$subcoho==1 & !is.na(S2$inz_mi))|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi")]
# S2.sub = na.omit(S2.sub)

require(survival)
rst = NULL;
##using Prentice or Lin Ying's weighting method
# for (m in c(S2_valid_measures)){
#   data.S2$metabolite = scale(log(data.S2[, m]))
#   model = cch(Surv(mi_time, inz_mi) ~ metabolite
#               + scale(ctalteru) #+ as.factor(ccsex)
#               + scale(ctbmi)## model 1
#               + as.factor(my.diab)  ##model 2
#               + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla) ##model 3+ total2HDL
#               + scale(log(cl_crp))  ##model 4
#               ,data = data.S2
#               , subcoh = ~subcoho
#               , id =~zz_nr
#               #, stratum= ~ as.factor(ccsex)
#               , method = "Prentice"
#               , cohort.size = 3244)
#   rst = rbind(rst, summary(model)$coefficients[1,])
# }
# #rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
# rownames(rst) = c(S2_valid_measures)
# #rst = cbind(rst, annotation[rownames(rst),])
# write.csv(rst, file = "metabolites_MI survival analysis_full model_plus medication_replication S2.csv")

## analysis using survey package
# require(survey)
# S2$Arg.Trp = S2$Arg/S2$Trp
# S2$eventrec = rep(0, nrow(S2))
# S2.cases = subset(S2, inz_mi==1)
# S2.cases$eventrec<-1
# S2.expd <- rbind(subset(S2, subcoho==1), S2.cases)
# S2.expd$stop <- with(S2.expd, ifelse(inz_mi&!eventrec, mi_time-1, mi_time))
# S2.expd$start <- with(S2.expd, ifelse(inz_mi&eventrec, mi_time-1, 0))
# S2.expd$event <- with(S2.expd, ifelse(inz_mi&eventrec, 1, 0))
# S2.expd$pwts[which(S2.expd$event==1)]=1
# #S2.expd$pwts[which(S2.expd$event==0&S2.expd$ccsex==1)]=1545/443
# #S2.expd$pwts[which(S2.expd$event==0&S2.expd$ccsex==2)]=1699/370
# S2.expd$pwts[which(S2.expd$event==0)]=3244/813
# S2.expd= subset(S2.expd, prev_mi==0)
# 
# dBarlow <- svydesign(id=~zz_nr+eventrec, starta=~subcoho+inz_mi, 
#                      data=S2.expd, weigts= S2.expd$pwts, probs = 1/S2.expd$pwts)
# 
# svycoxph(Surv(start, stop, event) ~ scale(log(Arg))
#          + scale(ctalteru) + as.factor(ccsex)
#          + scale(ctbmi)## model 1
#          + as.factor(my.diab)  ##model 2
#          + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla) ##model 3
#          + scale(log(cl_crp)),
#          design=dBarlow
# )


##Barlow's weighting
S2$Arg.Trp = S2$Arg/S2$Trp
S2.sub = S2[which(S2$subcoho==1|S2$inz_mi02==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time02", "inz_mi02","subcoho","prev_mi02","ctantihy", "ctmstati", "cl_ldla", "ctcigreg")]
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

##
model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ 
#                  scale(log(Arg)) + 
#                  scale(log(lysoPC_a_C17_0)) + 
#                  scale(log(Trp)) + 
#                  scale(log(PC_aa_C32_2))+
#                 scale(log(PC_aa_C36_3))+
#                 scale(log(lysoPC_a_C18_2)) + 
#                 scale(log(SM_C24_1))+                
              + scale(ctalteru) + as.factor(ccsex)
              + scale(ctbmi)## model 1
              + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
              + scale(log(cl_crp))  ##model 4
              #+ as.factor(ctmstati)
              + cluster(as.factor(zz_nr))
              ,data = S2.sub
              #,subset = S2.sub$ctmstati!=1
              ,weights = weight
              ,method = "breslow"
)
write.csv(cbind(summary(model)$coef, exp(confint(model))), file = "estimates of confounders in S2_model 4.csv.csv")


## association analysis on all metabolites
rst = NULL
for (m in c(S2_valid_measures,"Arg.Trp")){
  S2.sub$metabolite = scale(log(S2.sub[, m]))
  model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ metabolite  
              + scale(ctalteru) + as.factor(ccsex)## model 1
              + scale(ctbmi) + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
              + scale(log(cl_crp))  ##model 4
              #+ as.factor(ctmstati)
              #+ cluster(as.factor(zz_nr))
                ,data = S2.sub
                #,subset = S2.sub$ctmstati!=1
              ,weights = weight
                ,method = "breslow"
              )
  rst = rbind(rst, summary(model)$coefficients[1,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rst$lower = exp(rst$coef - 1.96*rst$se.coef.)
rst$upper = exp(rst$coef+1.96*rst$se.coef.)
rownames(rst) = c(S2_valid_measures,"Arg_Trp")
write.csv(rst, file = "metabolites_MI survival analysis_model 3_2002 S2 case cohort2_adj statines.csv")

require(gplots)
pdf("quintile plot of relative risk (ynorm)_model 4_S2.pdf", width = 12, height = 12)
par(mfrow =c(2,2));
yrange = NULL; RRquin = NULL
rst= NULL; rst2 = NULL
for(m in S2_valid_measures){
  m.conc=S2.sub[, m]
  metabo.quintile = cut(m.conc, breaks = quantile(m.conc, probs = seq(0, 1, 0.20), na.rm=T), include.lowest = T,ordered_result = F)
  model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ metabo.quintile +
                  scale(ctalteru) + as.factor(ccsex)
                + scale(ctbmi) + as.factor(my.diab)  ##model 2
                + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)
                + scale(log(cl_crp))
                ,data = S2.sub
                ,weights = weight
                ,method = "breslow"
                )
  rst = summary(model)$coefficients[1:4, ]
  rst2 = rbind(rst2, summary(model)$coefficients[1:4, ])

  #	interval = paste(round(exp(rst[,1] - 1.96* rst[,3]),3), 
  #			round(exp(rst[,1] + 1.96*rst[,3]), 3),
  #			sep = ",")
  #	interval = paste(round(rst[,2],3), interval, sep = "(")
  #	RRquin = cbind(RRquin, interval)
  upper = abs(exp(rst[,1] + rst[,3]) - rst[,2])
  lower = abs(exp(rst[,1] - rst[,3]) - rst[,2]) 
  if(max(rst[, 2]+upper)>1){
    yrange = c(min(rst[, 2]-upper), max(rst[, 2]+upper))
  }
  else{
    yrange = c( min(rst[, 2]-upper), 1)
  }
  x= tapply(m.conc, INDEX = metabo.quintile, median)
  plotCI(x = x[-1], y = rst[, 2], uiw = upper, liw = lower, main = m, 
         xlim = range(x), 
         ylim = yrange, 
         #xaxt = "n",
         log="y",
         pch=22,cex=3,pt.bg="black",
         #labels = levels(metabo.quintile), 
         ylab = "relative risk", xlab = "quintiles of metabolites (ratios)" )
  plotCI(x=x[1], y=1, uiw = 0, add=T, pch=21, cex=3, pt.bg="black")
  #axis(1, at = x, labels = levels(metabo.quintile), col.axis = "blue")
  lines(lowess(x, c(1, rst[,2]), f = 0.8), col = "red")
  abline(h = 1, lty = 2)
}
dev.off()
rownames(rst2) = rep(c(S2_valid_measures), each=4)
write.csv(rst2, file = "metabolites categorical_MI survival analysis_model 1_S2 case cohort.csv")


#test the trend
rst= NULL;
for(m in c(S2_valid_measures)){
  m.conc=scale(log(S2.sub[, m]))
  metabo.quintile = cut(m.conc, breaks = quantile(m.conc, probs = seq(0, 1, 0.20), na.rm=T), include.lowest = T,ordered_result = F)

  metabo.quintile.value = tapply(scale(log(S2.sub[, m])), INDEX = metabo.quintile, median, na.rm=T)		
  S2.sub.quintilevalue = rep(0, length(metabo.quintile))
  for (q in levels(metabo.quintile)){
    S2.sub.quintilevalue[which(metabo.quintile %in% q)]=metabo.quintile.value[q]
  }
  metabo.quintile = S2.sub.quintilevalue
  model = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ metabo.quintile +
                  scale(ctalteru) + as.factor(ccsex)
                + scale(ctbmi) + as.factor(my.diab)  ##model 2
                + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)
                + scale(log(cl_crp))
                #+cluster(as.factor(zz_nr))
                ,data = S2.sub
                ,weights = weight
                ,method = "breslow"
  )
  rst = rbind(rst, summary(model)$coefficients[1, ])
}
rownames(rst) = c(S2_valid_measures)
write.csv(rst, "metabolites categorical_test for trend_model 4_S2 case cohort.csv")

## By a matched case control study
# data = S2[which(S2$prev_mi==0),c("ctalteru", "ccsex", "ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola", "cl_hdla", "cl_ldla", "cl_crp", "inz_mi", "mi_time", S2_valid_measures)]
# data = na.omit(data)
# data$my.cigreg = as.numeric(data$my.cigreg)
# data$my.diab = as.numeric(data$my.diab)
# data$Arg_Trp = data$Arg/data$Trp
# 
# require(Matching)
# tr = data$inz_mi
# x = data[,c("ctalteru", "ccsex", "ctbmi", "cl_hdla","cl_chola")]# +ctbmi+cl_hdla+cl_chola
# glm1 <- glm(inz_mi~ctalteru+ccsex, family=binomial, data)
# rst.match = Match(Tr = tr, X =glm1$fitted.values, M =1, ties = T, replace=F, distance.tolerance = 0.5)
# 
# data = data[c(rst.match$index.control, unique(rst.match$index.treated)), ]
# data$ID = c(rst.match$index.treated,unique(rst.match$index.treated))
# 
# t.test(log(data$Arg[which(data$inz_mi==1)]), log(data$Arg[which(data$inz_mi==0)]), paired=T)
# 
# require(survival)
# rst = NULL
# for(m in S2_valid_measures){
#   data$m = scale(log(data[,m]))
#   model = clogit(inz_mi ~ m #+ ctalteru
#                  #+ scale(ctbmi)## model 1
#                  #+ scale(ctsysmm) #+ scale(cl_chola) + scale(cl_hdla) ##model 3
#                  #+ as.factor(my.diab)+ as.factor(my.cigreg)+as.factor(my.alkkon)
#                  #+ scale(cl_crp)  ##model 4
#                  + strata(ID)
#                  #, weights = data$weight
#                  ,data)
#   rst = rbind(rst, summary(model)$coef[1,])
# }
# rownames(rst) = S2_valid_measures
# 
# ### association with CRP
# rst = NULL
# for(i in S2_valid_measures){
#   S2$metabolite = scale(log(S2[,i]))
#   model = lm(scale(log(cl_crp)) ~ metabolite 
#              + scale(ctalteru) + as.factor(ccsex)
#              + scale(ctbmi)## model 1
#              + as.factor(my.diab)  ##model 2
#              + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)
#              , data = S2[which((S2$prev_mi==0&S2$subcoho==1)|S2$inz_mi==1),]
#              )
#   rst = rbind(rst, summary(model)$coef[2,])
# }
# rownames(rst) = S2_valid_measures
# write.csv(rst, "association with CRP_unadjsted_S2_case cohort.csv")
# 
# 
# model = lm(scale(log(cl_crp)) ~ scale(log(Arg)) + scale(log(Trp)) + scale(log(lysoPC_a_C17_0)) + scale(log(PC_aa_C32_2)) + scale(log(SM_C24_1)) + 
#            scale(ctalteru) + as.factor(ccsex)
#            + scale(ctbmi)## model 1
#            + as.factor(my.diab)  ##model 2
#            + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)
#            , data = S2[which((S2$prev_mi==0&S2$subcoho==1)|S2$inz_mi==1),]
# )

### association of metabolites with age
plot(log(Trp)~ctalteru, S2)
abline(lm(log(Trp)~ctalteru, S2))

model = lm(Trp~ctalteru + as.factor(ccsex) + scale(ctbmi) + scale(cl_hdla) + scale(cl_chola), S2)



