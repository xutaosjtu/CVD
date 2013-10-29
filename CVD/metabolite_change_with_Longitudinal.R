changeinBMI = F4[Cohort$zz_nr_f4,"utbmi"] - S4[Cohort$zz_nr_s4,"ltbmi"]
changeinLPC = log(F4[Cohort$zz_nr_f4,"C0"]) - log(S4[Cohort$zz_nr_s4,"C0"])


data$ltbmi.category = cut(data$ltbmi, breaks = c(0,25,60))

require(nlme)
rst=NULL
for(i in valid_measures){
  data$m = c(scale(log(S4[Cohort$zz_nr_s4, i])), scale(log(F4[Cohort$zz_nr_f4, i])))
  #data$m = scale(data$m)
  mixed.dum <- lme( m ~ as.factor(ltbmi.category)
                    + ltalteru + as.factor(lcsex) # crude model
                    + my.cigreg + my.alkkon 
                    + my.diab + ll_chola + ll_hdla +log(lh_crp) # multivariate model
                    ,random = ~  1 | participants, na.action=na.exclude, 
                    subset = which(data$ltantihy!=1),
                    data=data[order(data$participants),])
  rst = rbind(rst, summary(mixed.dum)$tTable[2,])
}
rst=data.frame(rst,
               fdr=p.adjust(rst[, 5], method="fdr"),
               bonf=p.adjust(rst[, 5], method="bonferroni")
)
rownames(rst)=valid_measures


mixed <- lm( scale(log(lysoPC_a_C18_2)) ~ ltbmi
                  + ltalteru + as.factor(lcsex) # crude model
                  + my.cigreg + my.alkkon 
                  + my.diab + ll_chola + ll_hdla +log(lh_crp) # multivariate model
                  , na.action=na.exclude, 
                  data = S4
                  )



#### longitudinal change of metabolites in MI patients
data$group = rep(interaction(data$my.mi[1:1002], data$my.mi[1003:2004]),2)

data$m = c(log(S4[Cohort$zz_nr_s4, i]), log(F4[Cohort$zz_nr_f4, i]))
# plotmeans(m~platform, data = data, subset = which(data$group=="(0,25].(0,25]"))
# plotmeans(m~platform, data = data, subset = which(data$group=="(25,60].(0,25]"), add = T)
# plotmeans(m~platform, data = data, subset = which(data$group=="(0,25].(25,60]"))
# plotmeans(m~platform, data = data, subset = which(data$group=="(25,60].(25,60]"), add = T)

pdf("lysoPCs change prospectively in the MI group.pdf")
par(mfrow = c(3,3))
for(i in valid_measures[c(20,30,33:40,45, 119)]){
#   S4$m = log(S4[,i])
#   model.S4 = lm(m ~ ltalteru + as.factor(lcsex) # crude model
#      + my.cigreg + my.alkkon + ltbmi
#      #+ my.diab + ll_chola + ll_hdla +log(lh_crp) # multivariate model
#      , data = S4)
#   S4$m[-na.action(model.S4)] = model.S4$residuals
#   S4$m[na.action(model.S4)] = NA
#   
#   F4$m = log(F4[,i])
#   model.F4 = lm(m ~ utalter + as.factor(ucsex) # crude model
#                 + my.cigreg + my.alkkon + utbmi
#                 #+ my.diab + ul_chola + ul_hdla +log(uh_crp) # multivariate model
#                 , data = F4)
#   F4$m[-na.action(model.F4)] = model.F4$residuals
#   F4$m[na.action(model.F4)] = NA
  
  data$m = c(S4[Cohort$zz_nr_s4, i], F4[Cohort$zz_nr_f4, i])
  
  model = lm(m~ ltalteru + as.factor(lcsex) # crude model
             + my.cigreg + my.alkkon + ltbmi
            # + my.diab + ll_chola + ll_hdla +log(lh_crp) # multivariate model
             , data = data)
  data$m[-na.action(model)] = model$residuals
  data$m[na.action(model)] = NA
  
  plotmeans(m ~ platform, data = data, subset = which(data$group == "0.1"), main = i)
}
dev.off()


require(nlme)
rst=NULL
for(i in valid_measures){
  data$m = c((log(S4[Cohort$zz_nr_s4, i])), (log(F4[Cohort$zz_nr_f4, i])))
  #data$m = scale(data$m)
  mixed.dum <- lme( m ~ as.factor(my.mi)
                    + ltalteru + as.factor(lcsex) # crude model
                    + my.cigreg + my.alkkon + ltbmi
                    + my.diab + ll_chola + ll_hdla +log(lh_crp) # multivariate model
                    ,random = ~  1 | participants, na.action=na.exclude, 
                    #subset = which(data$ltantihy!=1),
                    data=data[order(data$participants),])
  rst = rbind(rst, summary(mixed.dum)$tTable[2,])
}
rst=data.frame(rst,
               fdr=p.adjust(rst[, 5], method="fdr"),
               bonf=p.adjust(rst[, 5], method="bonferroni")
)
rownames(rst)=valid_measures
