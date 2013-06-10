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

###########################	Stroke	########################################
require(survival)
rst = NULL;
for (m in S4_valid_measures){
	metabolite = S4[, m]
	model = coxph(Surv(apo_time, inz_apo) ~  log(metabolite) +
					ltalteru + as.factor(lcsex) + ltbmi## model 1
			+(lp_diab_who06==4|lp_diab_who06==5)  ##model 2
			+log(ll_choln)+log(ll_hdln)+log(ltsysmm)+ as.factor(ltcigreg) + ltalkkon ##model 3
	 		#+ log(lh_crp) + total2HDL ##model 4
			#+ltdiamm ## model 5
			,subset = which(S4$prev_apo == 0&!(S4$apo_typ %in% c(5,1,2,9))),#
			S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = "Stroke survival analysis_model3.csv")

table(S4$inz_apo==1, S4$apo_typ)
S4$my.apo_typ = S4$apo_typ
S4$my.apo_typ[which(S4$apo_typ == 0)] = 0
S4$my.apo_typ[which(S4$apo_typ == 2)] = 1
S4$my.apo_typ[which(S4$apo_typ == 3)] = 2
S4$my.apo_typ[which(S4$apo_typ == 4)] = 2
S4$my.apo_typ[which(S4$apo_typ == 5)] = 3
S4$my.apo_typ[which(S4$apo_typ == 9)] = 4

require(caret)
require(pls)
require(gplots)
index = which((S4$my.apo_typ == 1|S4$my.apo_typ == 2)&S4$prev_apo==0)#
S4pls<-plsr(S4$my.apo_typ[index] ~ . , data=log2(S4[index, c(S4_valid_measures)]),  validation = "CV")
S4pls<-plsda(x=log2(S4[index, c(S4_valid_measures)]), y = as.factor(S4$my.apo_typ[index]))
#S4pls<-plsda(x=data.normalized, COPD.data$COPD, ncomp = 10)

color = greenred(24)[c(c(4:1)*2,c(18,20, 22, 24)) ]
plot(S4pls$scores[,c(1,2)],col = color[c(1,8)][S4$my.apo_typ[index]], pch=c(17, 19)[S4$my.apo_typ[index]])
legend(1, -3, 
		legend = c("Ischemic","Hemorrhagic"), 
		col = c(1:8), pch = c(17,19)
)

S4pca = prcomp(log(S4[index, S4_valid_measures]) )
S4pca = pcr(S4$my.apo_typ[index] ~ ., data = log2(S4[index, c(S4_valid_measures)]))
plot(S4pca$scores,col = color[S4$my.apo_typ[which(S4$my.apo_typ == 1|S4$my.apo_typ == 2)]], pch=c(17, 19)[S4$my.apo_typ[which(S4$my.apo_typ == 1|S4$my.apo_typ == 2)]])

difference = scan(what = character())
PC_aa_C28_1
PC_ae_C38_1
SM__OH__C22_1
SM_C16_0
SM_C24_0
C0
C10_2
Pro
Taurine
PC_ae_C40_4


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

require(survival)
rst = NULL; rst1 = NULL
rst2 = NULL; rst3 = NULL
for (m in S4_valid_measures){
	S4$metabolite = scale(log(S4[, m]))
	model = coxph(Surv(mi_time, inz_mi) ~ metabolite + #as.factor(ltnuecht) +
					scale(ltalteru) + as.factor(lcsex)
					+ scale(ltbmi)## model 1
					+ my.diab  ##model 2
					+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
	 				+ scale(lh_crp)  ##model 4
					#+ as.factor(ltmstati)& S4$ltmstati !=1
					,subset = which(S4$prev_mi==0),
					data = S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
	#rst1 = rbind(rst , summary(model)$coefficients[10,])
	#rst2 = rbind(rst2, summary(model)$coefficients[11,])
	#rst3 = rbind(rst3, summary(model)$coefficients[12,])
	#table(model$y[,2]) #number of sample used exactly in the estimation.
}
table(model$y[,2])
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = S4_valid_measures
#rst = cbind(rst, annotation[rownames(rst),])
write.csv(rst, file = "metabolites ratios_MI survival analysis_full model.csv")

plot(survfit(Surv(mi_time, S4$inz_mi)~(log(S4$PC_aa_C32_2) > 1.2), S4, subset= which(S4$prev_mi == 0)), log = "y", col = c("red","green"))

############	Hazardous ratios in different quantiles	################
require(gplots)
pdf("quintile plot of replative risk (ynorm)_full model _decile.pdf", width = 12, height = 12)
par(mfrow =c(2,2));
yrange = NULL; RRquin = NULL
for(m in candidates[1:16]){
	m.conc=S4[, m]
	metabo.quintile = cut(m.conc, breaks = 
					#range(exp(S4[, m]))[1] + abs(range(exp(S4[, m]))[1]-range(exp(S4[, m]))[2])/6*(1:6), 
					quantile(m.conc, probs = seq(0, 1, 0.1)), 
			include.lowest = T,ordered_result = F)
	model1 = coxph(Surv(mi_time, inz_mi) ~ metabo.quintile +
					scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
					#+ my.diab  ##model 2
					#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
					#+ scale(lh_crp)  ##model 4&S4$ltmstati !=1
			,subset = which(S4$prev_mi == 0),
			S4)
	rst = summary(model1)$coefficients[1:9, ]
	
	model2 = coxph(Surv(mi_time, inz_mi) ~ metabo.quintile +
					scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
			+ my.diab  ##model 2
			+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3
			+ scale(lh_crp)  ##model 4&S4$ltmstati !=1
			,subset = which(S4$prev_mi == 0),
			S4)
	rst2=summary(model2)$coefficients[1:9,]
	
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

#test the trend
rst= NULL;
for(m in metabo.selected3){
	metabo.quintile = cut(S4[, m], breaks = quantile(S4[, m], probs = seq(0, 1, 0.2)), include.lowest = T,ordered_result = F)
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
					scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
					#+ my.diab  ##model 2
					#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
					#+ scale(lh_crp)  ##model 4
			,subset = which(S4$prev_mi == 0&S4$ltmstati !=1),
			S4)
	rst = rbind(rst, summary(model)$coefficients[1, ])
}


###metabolites associated with MI
metabo.asso = scan(what = character())
Arg
Trp
lysoPC_a_C16_0
lysoPC_a_C17_0
lysoPC_a_C18_1
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2
PC_aa_C34_2
PC_aa_C34_3
PC_aa_C36_2
PC_aa_C36_3
PC_ae_C36_2
PC_ae_C38_0
PC_ae_C38_2
PC_ae_C40_1



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



