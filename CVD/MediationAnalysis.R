S4$metabolite = scale(log(S4[, m]))
model.1 = coxph(Surv(mi_time, inz_mi)  ~  scale((lh_crp))##model 4
              + scale(ltalteru) + as.factor(lcsex)## model 1
              + scale(ltbmi)+ as.factor(my.diab)  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
              ,subset = which(S4$prev_mi==0&!is.na(S4$inz_mi)),#
              data = S4)

model.2.coef = NULL
for(m in c("Arg", "lysoPC_a_C17_0", "lysoPC_a_C18_2")){
  S4$m = scale(log(S4[,m]))
  model.2 = lm(m ~ scale((lh_crp))##model 4
               + scale(ltalteru) + as.factor(lcsex)## model 1
               + scale(ltbmi)+ as.factor(my.diab)  ##model 2
               + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
               ,subset = which(S4$prev_mi==0&!is.na(S4$inz_mi)),#
               data = S4)
  model.2.coef = rbind(model.2.coef, summary(model.2)$coef[2,])
}
rownames(model.2.coef) = c("Arg", "lysoPC_a_C17_0", "lysoPC_a_C18_2")

require(survival)
model.3 = coxph(Surv(mi_time, inz_mi) ~ scale(log(S4[,"Arg"])) + scale(log(S4[,"lysoPC_a_C17_0"])) + scale(log(S4[,"lysoPC_a_C18_2"])) 
              + scale((lh_crp))
              + scale(ltalteru) + as.factor(lcsex)## model 1
              + scale(ltbmi)+ as.factor(my.diab)  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
              ,subset = which(S4$prev_mi==0),#
              data = S4)



## S2
S2.sub$metabolite = scale(log(S2.sub[, m]))
model.1 = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ 
                scale((cl_crp))  
              + scale(ctalteru) + as.factor(ccsex)## model 1
              + scale(ctbmi) + as.factor(my.diab)  ##model 2
              + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
                ,data = S2.sub
                ,weights = weight
                ,method = "breslow"
)

model.2.coef = NULL
for(m in c("Arg", "lysoPC_a_C17_0", "lysoPC_a_C18_2")){
  S2$m = scale(log(S2[,m]))
  model.2 = lm(m ~ scale((cl_crp))  
                + scale(ctalteru) + as.factor(ccsex)## model 1
                + scale(ctbmi) + as.factor(my.diab)  ##model 2
                + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
                ,data = S2
  )
  model.2.coef = rbind(model.2.coef, summary(model.2)$coef[2,])
}
rownames(model.2.coef) = c("Arg", "lysoPC_a_C17_0", "lysoPC_a_C18_2")

require(survival)
model.3 = coxph(Surv(mi_time.start, mi_time.end, inz_mi02) ~ 
                  scale(log(S2.sub[,"Arg"])) + scale(log(S2.sub[,"lysoPC_a_C17_0"])) + scale(log(S2.sub[,"lysoPC_a_C18_2"])) 
                + scale((cl_crp))  
                + scale(ctalteru) + as.factor(ccsex)## model 1
                + scale(ctbmi) + as.factor(my.diab)  ##model 2
                + scale(ctsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(cl_chola) + scale(cl_hdla) ##model 3
                
                ,data = S2.sub
                ,weights = weight
                ,method = "breslow")
