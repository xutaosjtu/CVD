S4$metabolite = scale(log(S4[, m]))
model = lm(metabolite ~  scale(log(lh_crp))##model 4
              + scale(ltalteru) + as.factor(lcsex)## model 1
              + scale(ltbmi)+ as.factor(my.diab)  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
              ,subset = which(S4$prev_mi==0&!is.na(S4$inz_mi)),#
              data = S4)

S4$metabolite = scale(log(S4[, m]))
model = survreg(Surv(mi_time, inz_mi) ~ metabolite 
              + scale(ltalteru) + as.factor(lcsex)## model 1
              + scale(ltbmi)+ as.factor(my.diab)  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
              + scale(log(lh_crp))##model 4
              #+ as.factor(ltmstati)
              #+ as.factor(ltnuecht)
              ,subset = which(S4$prev_mi==0),#
              data = S4)



