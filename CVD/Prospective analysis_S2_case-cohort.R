S2$Arg.Trp = S2$Arg/S2$Trp
tmp = S2[which((S2$subcoho==1 & !is.na(S2$inz_mi))|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi")]
tmp = na.omit(tmp)

require(survival)
rst = NULL;
##using Prentice or Lin Ying's weighting method
for (m in c(S2_valid_measures,"Arg.Trp")){
  tmp$metabolite = scale(log(tmp[, m]))
  model = cch(Surv(mi_time, inz_mi) ~ metabolite + #as.factor(ltnuecht) +
                scale(ctalteru) #+ as.factor(ccsex)
              + scale(ctbmi)## model 1
              + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla) ##model 3+ total2HDL
              + scale(log(cl_crp))  ##model 4
              ,data = tmp
              , subcoh = tmp$subcoho
              , id =~zz_nr
              #, stratum= ~ as.factor(ccsex)
              , method = "Prentice"
              , cohort.size = 3408)
  rst = rbind(rst, summary(model)$coefficients[1,])
  #rst1 = rbind(rst , summary(model)$coefficients[10,])
  #rst2 = rbind(rst2, summary(model)$coefficients[11,])
  #rst3 = rbind(rst3, summary(model)$coefficients[12,])
  #table(model$y[,2]) #number of sample used exactly in the estimation.
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = c(S2_valid_measures,"Arg_Trp")
#rst = cbind(rst, annotation[rownames(rst),])
write.csv(rst, file = "metabolites_MI survival analysis_full model_plus medication_replication S2.csv")

##Barlow's weighting
S2$Arg.Trp = S2$Arg/S2$Trp
tmp = S2[which(S2$subcoho==1|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho","prev_mi","ctantihy")]
#tmp = na.omit(tmp)
tmp = tmp[which(tmp$prev_mi!=1),]
tmp$mi_time.start = 0
tmp$mi_time.start[which(tmp$inz_mi==1)] = tmp$mi_time[which(tmp$inz_mi==1)]-1
tmp$mi_time.end = tmp$mi_time

#cases in subcohort
#tmp$atcontrol = tmp$inz_mi
subcohort.cases = tmp[which(tmp$subcoho==1&tmp$inz_mi==1),]
#subcohort.cases$atcontrol = 0
subcohort.cases$inz_mi=0
subcohort.cases$mi_time.start = 0
subcohort.cases$mi_time.end = subcohort.cases$mi_time.end-1

tmp = rbind(tmp, subcohort.cases)

weight=rep(1,nrow(tmp))
weight[which(tmp$ccsex==1&tmp$subcoho==1 & tmp$inz_mi==0)]= 1679/306
weight[which(tmp$ccsex==2&tmp$subcoho==1 & tmp$inz_mi==0)]= 1729/289
tmp$weight = weight

rst = NULL
for (m in c(S2_valid_measures,"Arg.Trp")){
  tmp$metabolite = scale(log(tmp[, m]))
  model = coxph(Surv(mi_time.start, mi_time.end, inz_mi) ~ metabolite + 
                scale(ctalteru) + as.factor(ccsex)
              + scale(ctbmi)## model 1
              + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla) ##model 3
              + scale(log(cl_crp))  ##model 
              + cluster(zz_nr)
                ,data = tmp
              ,weights = weight
                ,ties = "breslow"
              )
  rst = rbind(rst, summary(model)$coefficients[1,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = c(S2_valid_measures,"Arg_Trp")
write.csv(rst, file = "metabolites_MI survival analysis_full model_S2 case cohort.csv")

require(gplots)
pdf("quintile plot of replative risk (ynorm)_full model _decile.pdf", width = 12, height = 12)
par(mfrow =c(2,2));
yrange = NULL; RRquin = NULL
for(m in candidates[1:16]){
  m.conc=tmp[, m]
  metabo.quintile = cut(m.conc, breaks = quantile(m.conc, probs = seq(0, 1, 0.2)), 
                        include.lowest = T,ordered_result = F)
  model = coxph(Surv(mi_time, inz_mi) ~ metabo.quintile +
                  scale(ctalteru) + as.factor(ccsex)
                + scale(ctbmi)## model 1
                + as.factor(my.diab)  ##model 2
                + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)+ scale(log(cl_crp))
                ,data = tmp
                ,weights = weight
                ,ties = "breslow"
                )
  rst = summary(model1)$coefficients[1:9, ]
  
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
  plotCI(x = x[2:10], y = rst[, 2], uiw = upper, liw = lower, main = m, 
         xlim = range(x), 
         ylim = yrange, 
         #xaxt = "n",
         #log="y",
         pch=22,cex=3,pt.bg="black",
         #labels = levels(metabo.quintile), 
         ylab = "relative risk", xlab = "quintiles of metabolites (ratios)" )
  plotCI(x=x[1], y=1, uiw = 0, add=T, pch=22, cex=3, pt.bg="black")
  plotCI(x = x[2:10], y = rst2[, 2], uiw = upper, liw = lower, main = m, 
         xlim = range(x), 
         ylim = yrange, 
         #xaxt = "n",
         #log="y",
         pch=21,cex=3,pt.bg="grey",
  )
  plotCI(x=x[1], y=1, uiw = 0, add=T, pch=21, cex=3, pt.bg="grey")
  #axis(1, at = x, labels = levels(metabo.quintile), col.axis = "blue")
  lines(lowess(x, c(1, rst[,2]), f = 0.8), col = "red")
  abline(h = 1, lty = 2)
}
dev.off()



## analysis using survey package
S2$Arg.Trp = S2$Arg/S2$Trp
S2$eventrec = rep(0, nrow(S2))
S2.cases = subset(S2, inz_mi==1)
S2.cases$eventrec<-1
S2.expd <- rbind(subset(S2, subcoho==1), S2.cases)
S2.expd$stop <- with(S2.expd, ifelse(inz_mi&!eventrec, mi_time-1, mi_time))
S2.expd$start <- with(S2.expd, ifelse(inz_mi&eventrec, mi_time-1, 0))
S2.expd$event <- with(S2.expd, ifelse(inz_mi&eventrec, 1, 0))
S2.expd$pwts[which(S2.expd$event==1)]=1
S2.expd$pwts[which(S2.expd$event==0&S2.expd$ccsex==1)]=1679/443
S2.expd$pwts[which(S2.expd$event==0&S2.expd$ccsex==2)]=1729/370
S2.expd= subset(S2.expd, prev_mi==0)

dBarlow <- svydesign(id=~zz_nr+eventrec, starta=~subcoho+inz_mi, 
                     data=S2.expd, weigts= S2.expd$pwts)

svycoxph(Surv(start, stop, event) ~ scale(log(PC_ae_C40_1))
         + scale(ctalteru) + as.factor(ccsex)
         + scale(ctbmi)## model 1
         + as.factor(my.diab)  ##model 2
         + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla) ##model 3
         + scale(log(cl_crp)),
         design=dBarlow
         )

