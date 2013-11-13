# The main focuse of the prospective analysis rest in two parts:
# first, is it possible to find better predictors other than the conventional biomarkers of CVD
# second, to see the propective changes of metabolite concentration in CVD patients.
# 
# Author: tao.xu
###############################################################################

rownames(S4) = S4$ZZ_nr.x
tmp = S4[which(S4[as.character(Cohort$zz_nr_s4),"ltjnc7"]>2),c("ltmbbl","ltmace","ltmata")]
medication = apply(tmp, 1, function(x) sum(x==1))
table(S4$ltmbbl[which(S4[,diseases[1]]>2)])
table(S4$ltmace[which(S4[,diseases[1]]>2)])
table(S4$ltmata[which(S4[,diseases[1]]>2)])

rownames(F4) = F4$ZZ_nr
tmp = F4[which(S4[as.character(Cohort$zz_nr_s4),"ltjnc7"]>2),c("utmbbl","utmace","utmata")]
medication = apply(tmp, 1, function(x) sum(x==1))

medication<-function(index){
	tmp = S4[index, c("ltmbbl","ltmace","ltmata")]
	m = apply(tmp, 1, function(x) sum(x==1))
	print("In S4, ")
	print(table(m))
	
	tmp = F4[index, c("utmbbl","utmace","utmata")]
	m = apply(tmp, 1, function(x) sum (x==1))
	print("In F4,")
	table(m)
}

index = which(S4[as.character(Cohort$zz_nr_s4),"ltjnc7"]<2&F4[as.character(Cohort$zz_nr_f4),"utjnc7"]>2)



############	calculate the residues	########
data = data.frame(log(S4[,valid_measures]),
		"alteru" = S4$ltalteru, 
		"sex" = as.factor(S4$lcsex),
		"bmi" = S4$ltbmi,
		"cigreg" = as.factor(S4$ltcigreg),
		"alkkon" = as.factor(S4$my.alkkon),
		"diabetes" = as.factor((S4$lp_diab_who06==4|S4$lp_diab_who06==5)),
		"chol" = log(S4$ll_chola),
		"HDL" = log(S4$ll_hdla),
		"med1" = as.factor(S4$ltmbbl),
		"med2" = as.factor(S4$ltmace),
		"med3" = as.factor(S4$ltmata),
		"hyper" = S4$lthyact
)
data = data[Cohort$zz_nr_s4, ]
S4.residue = residue(data, valid_measures, adj = colnames(data)[122:132], control_group = 1:dim(data)[1])

data = data.frame(log(F4[,valid_measures]),
		"alteru" = F4$utalteru, 
		"sex" = as.factor(F4$ucsex),
		"bmi" = F4$utbmi,
		"cigreg" = as.factor(F4$utcigreg),
		"alkkon" = as.factor(F4$my.alkkon),
		"diabetes" = as.factor((F4$uk_diab_who06==4|F4$uk_diab_who06==5)),
		"chol" = log(F4$ul_chola),
		"HDL" = log(F4$ul_hdla),
		"med1" = as.factor(F4$utmbbl),
		"med2" = as.factor(F4$utmace),
		"med3" = as.factor(F4$utmata),
		"hyper" = F4$uthyact
)
data = data[Cohort$zz_nr_f4, ]
F4.residue = residue(data, valid_measures, adj = colnames(data)[122:132], control_group = 1:dim(data)[1])

metabo.selected = scan(what = character())
C8
C10
C10_1
C14_2
Phe
Pro
lysoPC_a_C17_0
lysoPC_a_C18_2
PC_ae_C32_2
PC_ae_C34_1
PC_ae_C34_2
PC_ae_C34_3
PC_ae_C36_1
PC_ae_C36_2
PC_ae_C42_3
SM__OH__C16_1
SM__OH__C22_2
H1

pdf(file = "S4 F4 prospective change.pdf", width = 15, height = 27)
close.screen(all = T)
split.screen(c(6,3))
k = 1
for(i in metabo.selected){
	screen(k)
	data = data.frame("m" = c(S4.residue[,i], F4.residue[,i]),  
			"group" = rep(interaction(as.factor(S4[Cohort$zz_nr_s4, "lthyact"]),as.factor(F4[Cohort$zz_nr_f4, "uthyact"])),2),
			"t" = c(rep(1, dim(S4.residue)[1]),rep(2, dim(S4.residue)[1]))
	)
	plotmeans(m ~ t, data, subset = which(data$group ==1.1 ), ylim = c(-0.1,0.1), main = i)
	plotmeans(m ~ t, data, subset = which(data$group ==1.2 ), col = "blue", add = T)
	plotmeans(m ~ t, data, subset = which(data$group ==2.1 ), col = "red", add = T)
	plotmeans(m ~ t, data, subset = which(data$group ==2.2 ), col = "green", add = T)
	k = k+1
}
dev.off()


###################	Myocardial infarction ##############################
S4$my.cigreg = S4$ltcigreg
S4$my.cigreg[which(S4$ltcigreg ==2)] = 1
S4$my.cigreg[which(S4$ltcigreg ==3)] = 2
S4$my.cigreg[which(S4$ltcigreg ==4)] = 3
S4$my.cigreg = 3-S4$my.cigreg
S4$my.cigreg = factor(S4$my.cigreg, ordered = F)

S4$total2HDL = S4$ll_chola/S4$ll_hdla

plot(
		S4$mi_time[which(S4$prev_mi==0&!is.na(S4$inz_mi))],
		S4_raw$mi_time[which(!is.na(S4_raw$inz_mi)&S4_raw$ltnuecht==1)]
)

##
##  estimates of the confounders
model = coxph(Surv(mi_time, inz_mi) ~ 
                scale(log(Arg)) + 
                 scale(log(lysoPC_a_C17_0)) + 
                 scale(log(Trp)) + 
                 scale(log(PC_aa_C32_2))+
#                 scale(log(PC_aa_C36_3))+
                 scale(log(lysoPC_a_C18_2)) + 
#                 scale(log(SM_C24_1))+
                scale(ltalteru) + as.factor(lcsex)
              + scale(ltbmi) + as.factor(my.diab)  ##model 2
              + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3
              + scale(log(lh_crp))  ##model 4
              ,subset = which(S4$prev_mi==0),
              data = S4)
write.csv(cbind(summary(model)$coef, exp(confint(model))), file = "estimates of confounders plus original four metabolites in S4_model 4_5 metabolites.csv")
##

## association analysis
require(survival)
rst = NULL; rst1 = NULL
rst2 = NULL; rst3 = NULL
for (m in S4_valid_measures){
	S4$metabolite = scale(log(S4[, m]))
	model = coxph(Surv(mi_time, inz_mi) ~ metabolite 
					+ scale(ltalteru) + as.factor(lcsex)
					+ scale(ltbmi)## model 1
					+ as.factor(my.diab)  ##model 2
					+ scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
	 				+ scale(log(lh_crp))##model 4
	        #+ as.factor(ltmstati)
          #+ as.factor(ltantihy)
					,subset = which(S4$prev_mi==0 & S4$ltmstati !=1),#
					data = S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
	#rst1 = rbind(rst , summary(model)$coefficients[2,])
	#rst2 = rbind(rst2, summary(model)$coefficients[14,])
	#rst3 = rbind(rst3, summary(model)$coefficients[12,])
	#table(model$y[,2]) #number of sample used exactly in the estimation.
}
table(model$y[,2])
rst = data.frame(rst,FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rst$lower = exp(rst$coef - 1.96*rst$se.coef.)
rst$upper = exp(rst$coef+1.96*rst$se.coef.)
rownames(rst) = S4_valid_measures
#rst = cbind(rst, annotation[rownames(rst),])
write.csv(rst, file = "metabolites_MI survival analysis_S4_unadj.csv")
rst2 = data.frame(rst2, FDR = p.adjust(rst2[,5], method = "BH"), bonferroni = p.adjust(rst2[,5], method = "bonferroni"))
rownames(rst2) = S4_valid_measures
write.csv(rst2, file = "metabolites fasting interaction_MI survival analysis_S4_full model_with nonfasting.csv")

index=sapply(S2[,S2_valid_measures], function(x) which(abs(x)>mean(x,na.rm=T)+5*sd(x,na.rm=T)|abs(x)<mean(x,na.rm=T)-5*sd(x,na.rm=T))) 
for(i in S2_valid_measures){
  S2[index[[i]],i]=NA
}

S2$Arg.Trp = S2$Arg/S2$Trp
tmp = S2[which((S2$subcoho==1 & !is.na(S2$inz_mi))|S2$inz_mi==1),c("ctalteru", "ccsex","ctbmi","my.diab","ctsysmm","my.cigreg","my.alkkon","cl_chola","cl_hdla","cl_crp",S2_valid_measures,"Arg.Trp","zz_nr","mi_time", "inz_mi","subcoho")]
tmp = na.omit(tmp)

plot(survfit(Surv(mi_time, S4$inz_mi)~(log(S4$PC_aa_C32_2) > 1.2), S4, subset= which(S4$prev_mi == 0)), log = "y", col = c("red","green"))

############	Hazardous ratios in different quantiles	################
require(gplots)
pdf("quartile plot of replative risk (ynorm)_model 1_S4.pdf", width = 12, height = 12)
par(mfrow =c(2,2));
yrange = NULL; RRquin = NULL
rst = NULL; rst2 = NULL
for(m in S4_valid_measures){
	m.conc=S4[, m]
	metabo.quintile = cut(m.conc, breaks = quantile(m.conc, probs = seq(0, 1, 0.25), na.rm = T), 
			include.lowest = T,ordered_result = F)
	model1 = coxph(Surv(mi_time, inz_mi) ~ metabo.quintile +
					scale(ltalteru) + as.factor(lcsex) 
         #+ scale(ltbmi) + my.diab  ##model 2
				 #+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
				 #+ scale(log(lh_crp))  ##model 4&S4$ltmstati !=1
         
			,subset = which(S4$prev_mi == 0),
			S4)
	rst = summary(model1)$coefficients[1:3, ]
  rst2 = rbind(rst2, summary(model1)$coefficients[1:3, ])
  
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
	plotCI(x = x[2:length(x)], y = rst[, 2], uiw = upper, liw = lower, main = m, 
			xlim = range(x), 
			ylim = yrange, 
			#xaxt = "n",
			#log="y",
			pch=22,cex=3,pt.bg="black",
			#labels = levels(metabo.quintile), 
			ylab = "relative risk", xlab = "quintiles of metabolites (ratios)" )
	plotCI(x=x[1], y=1, uiw = 0, add=T, pch=22, cex=3, pt.bg="black")
  
	lines(lowess(x, c(1, rst[,2]), f = 0.8), col = "red")
	abline(h = 1, lty = 2)
}
dev.off()
rownames(rst2) = rep(c(S4_valid_measures), each=4)
write.csv(rst2, file = "metabolites categorical_MI survival analysis_model 1_S4.csv")




#test the trend
rst= NULL;
for(m in S4_valid_measures){
	metabo.quintile = cut(S4[, m], breaks = quantile(S4[, m], probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T,ordered_result = F)
	#log(S4[,m])
	if(m == "Arg.Trp"){
		metabo.quintile.value = tapply(scale(S4[, m]), INDEX = metabo.quintile, median, na.rm=T)
	}
	else{
		metabo.quintile.value = tapply(scale(log(S4[, m])), INDEX = metabo.quintile, median, na.rm=T)		
	}
	tmp = rep(0, length(metabo.quintile))
	for (q in levels(metabo.quintile)){
		tmp[which(metabo.quintile %in% q)]=metabo.quintile.value[q]
	}
	metabo.quintile = tmp
	model = coxph(Surv(mi_time, inz_mi) ~ metabo.quintile +
					scale(ltalteru) + as.factor(lcsex) 
          #+ scale(ltbmi) + my.diab  ##model 2
					#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
					#+ scale(log(lh_crp))  ##model 4
			,subset = which(S4$prev_mi == 0),
			S4)
	rst = rbind(rst, summary(model)$coefficients[1, ])
}
rownames(rst) = S4_valid_measures
write.csv(rst, "metabolites categorical_test for trend_model 1_S4.csv")

###metabolites associated with MI
metabo.asso = scan(what = character())
Arg
Trp
lysoPC_a_C16_0
lysoPC_a_C16_1
lysoPC_a_C17_0
lysoPC_a_C18_0
lysoPC_a_C18_1
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2
PC_aa_C34_2
PC_aa_C34_3
PC_aa_C36_2
PC_aa_C36_3
PC_ae_C36_1
PC_ae_C36_2
PC_ae_C38_0
PC_ae_C40_1
SM_C24_1




metabo.ratio.asso = scan(what = character())
Arg.Trp
Arg.lysoPC_a_C16_0
Arg.lysoPC_a_C17_0
Arg.lysoPC_a_C18_2
Arg.PC_aa_C28_1
Arg.PC_aa_C32_2
Arg.PC_aa_C32_3
Arg.PC_aa_C34_2
Arg.PC_aa_C34_3
Arg.PC_aa_C36_2
Arg.PC_aa_C36_3
Arg.PC_ae_C36_1
Arg.PC_ae_C36_2
Arg.PC_ae_C38_2
Arg.PC_ae_C40_1
PC_aa_C32_2.PC_aa_C32_3
PC_aa_C32_2.PC_aa_C34_2
PC_aa_C32_2.PC_aa_C34_3
PC_aa_C32_2.PC_aa_C36_2
PC_aa_C32_2.PC_aa_C36_3





clinical = c("ltalteru", "ltbmi", "lcsex","my.diab", "ltsysmm", "ll_hdln", "ll_choln", "my.cigreg", "my.alkkon", "lh_crp")



