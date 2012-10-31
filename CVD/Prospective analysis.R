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


###	new developped Myocardial Infarction
tmps4 = as.factor(S4[Cohort$zz_nr_s4, diseases.s4[i]])
tmpf4 = as.factor(F4[Cohort$zz_nr_f4, diseases.f4[i]])
f = interaction(tmpf4, tmps4)

CVD.cox = coxph(Surv(t_time, CVD) ~  log(metabolite) + ltalteru + log(ltsysmm) + log(ll_hdln) + log(ll_choln) + as.factor(lp_diab_who06) + as.factor(lcsex) , subset = which(f %in% c("2.1", "2.2")), data)


#index = which(F4[as.character(Cohort$zz_nr_f4),"CVD"])
#data.newCVD = data.frame(S4[as.character(Cohort$zz_nr_s4)[index], c("CVD", "ltalteru") ] , F4[as.character(Cohort$zz_nr_f4)[index],c("CVD","utmialt","utschalt","utalter")])
#write.csv(data.newCVD, file = "newly developped CVDs.csv")
#
#index2 = which(S4[as.character(Cohort$zz_nr_s4),"CVD"])
#data.newCVD = data.frame(S4[as.character(Cohort$zz_nr_s4)[index2], c("CVD", "ltalteru") ] , F4[as.character(Cohort$zz_nr_f4)[index2],c("CVD","utmialt","utschalt","utalter")])
#write.csv(data.newCVD, file = "CVDs in S4.csv")
#
#table(F4[as.character(Cohort$zz_nr_f4)[-c(index, index2)],"CVD"])
#table(S4[as.character(Cohort$zz_nr_s4)[-c(index, index2)],"CVD"])
#
## assessing the biomarkers for prediction of CVD (on 26 cases!!!) 
#index = which(!(S4[as.character(Cohort$zz_nr_s4),"CVD"]))
#substr(feature.cont, 1, 2) <- "l"
#substr(feature.disc, 1, 2) <- "l"
#rst = NULL;
#for (m in S4_valid_measures){
#	data = data.frame(
#			CVD = F4[as.character(Cohort$zz_nr_f4), "CVD"],
#			
#			CVDS4 = S4[as.character(Cohort$zz_nr_s4), "CVD"], 
#			
#			S4[as.character(Cohort$zz_nr_s4), m],
#			
#			S4[as.character(Cohort$zz_nr_s4),feature.cont],
#			
#			S4[as.character(Cohort$zz_nr_s4),feature.disc], 
#			
#			t_time = apply(F4[as.character(Cohort$zz_nr_f4), c("utmialt","utschalt", "utalteru")], 1, min, na.rm = T)-S4[as.character(Cohort$zz_nr_s4),"ltalteru"]
#	)
#	
#	data = na.omit(data)
	
#	index = which(!data$CVDS4 & data$t_time>0)
#	
#	colnames(data) = c("CVD" , "CVDS4" , "metabolite", "ltalteru" , "ltsysmm" , "ll_hdln" , "ll_choln" , "lp_diab_who06" , "lcsex" , "t_time")
#	
#	CVD.cox = coxph(Surv(t_time, CVD) ~  log(metabolite) + ltalteru + log(ltsysmm) + log(ll_hdln) + log(ll_choln) + as.factor(lp_diab_who06) + as.factor(lcsex) , subset = index , data)
#	
#	rst = rbind(rst , summary(CVD.cox)$coefficients[1,])
#}
#rownames(rst) = S4_valid_measures


###############################separate CVDs#######################################
# model: age, sex, bmi, alcohol consumption
#
#
diseases.s4 = scan(what = character())
ltjnc7
ltschl
ltmi
lc044f_1

#names(diseases) = c("Hypertension","Stroke","Myocardio","Heart failure")
#
#
diseases.f4 = scan(what = character())
utjnc7
utschl
utmi
us_c04a

#names(diseases) = c("Hypertension","Stroke","Myocardio","Heart failure")
#

table(S4[Cohort$zz_nr_s4, "CVD"], F4[Cohort$zz_nr_f4, "CVD"])
table(tmps4 = S4[Cohort$zz_nr_s4, "ltjnc7"], tmpf4 = F4[Cohort$zz_nr_f4, "utjnc7"])
table(tmps4 = S4[Cohort$zz_nr_s4, "ltschl"], tmpf4 = F4[Cohort$zz_nr_f4, "utschl"], useNA = "always")
table(tmps4 = S4[Cohort$zz_nr_s4, "ltmi"], tmpf4 = F4[Cohort$zz_nr_f4, "utmi"])
table(tmps4 = S4[Cohort$zz_nr_s4, "lc044f_1"], tmpf4 = F4[Cohort$zz_nr_f4, "us_c04a"])



for(i in 2:4){
	
	tmps4 = S4[Cohort$zz_nr_s4, diseases.s4[i]]
	tmpf4 = F4[Cohort$zz_nr_f4, diseases.f4[i]]
	feature = data.frame(tmps4, tmpf4)
	subset = intersect (which(feature[,1] == 2), which(feature[,2] !=3 ))
	rst = Comparison.prospective(S4[ Cohort$zz_nr_s4, ], feature, metabo = S4_valid_measures, subset = subset, adj = model3[[i]])
	write.csv(rst, file = paste(names(diseases.f4)[i], "prospective at S4 baseline_model3.csv"))
}

#	tmps4 = S4[Cohort$zz_nr_s4, diseases.s4[i]]
#	tmpf4 = F4[Cohort$zz_nr_f4, diseases.f4[i]]
##	pheno = interaction(tmps4, tmpf4)
##	pheno = which(pheno == )
#	rst = NULL
#	for(j in 1:length(S4_valid_measures)){		
#		data = data.frame(log(S4[Cohort$zz_nr_s4,S4_valid_measures[j]]), 
#						S4[Cohort$zz_nr_s4, c("ltalter", "lcsex", "ltbmi" , "ltalkkon")], tmps4, tmpf4 )
#		data = data[which(data$tmps4 == 2),]
#		data = data[which(data$tmpf4 != 3),]
##		data = data[which(data$tmps4),]
#		model = glm(interaction(tmps4, tmpf4) ~ ., data, family = binomial(link = "logit"))
#		rst = rbind(rst, summary(model)$coefficients[2,])
#	}
#	rownames(rst) = S4_valid_measures
#	
#	rst = apply(
#			S4[Cohort$zz_nr_s4, S4_valid_measures],
#			2,
#			function(x) tapply(x, INDEX = interaction(tmps4, tmpf4), mean)
#	)


###########################	Stroke	########################################
#require(nlme)
#
#valid_measures = intersect(S4_valid_measures, F4_valid_measures)
#participants=rep(1:1009,2)
#
#data=data.frame(
#		participants, 
#		disease =  as.factor(c(S4[Cohort$zz_nr_s4, "ltjnc7"], F4[Cohort$zz_nr_f4, "utjnc7"])),
#		rbind(as.matrix(S4[ Cohort$zz_nr_s4, feature.cont.s4]), as.matrix(F4[ Cohort$zz_nr_f4, feature.cont.f4])),
#		apply(rbind(as.matrix(S4[ Cohort$zz_nr_s4, feature.disc.s4]), as.matrix(F4[ Cohort$zz_nr_f4, feature.disc.f4])), 2, function(x) as.factor(x))
#)
#
#rst=NULL
#for(i in valid_measures){
#
#	m = c(S4[Cohort$zz_nr_s4, i], F4[Cohort$zz_nr_f4, i])
#	
#	mixed.dum <- lme( m ~  ltsysmm + ltbmi + ltalteru + ltalkkon + lttumf + waist2hip + ltrauchp + lcsex , random = ~  1 | participants, na.action=na.exclude, data=data)
#	
#	rst = rbind(rst, summary(mixed.dum)$tTable[5,])
#}
#rst=data.frame(rst,
#		fdr=p.adjust(rst[, 5], method="fdr"),
#		bonf=p.adjust(rst[, 5], method="bonferroni")
#)
#rownames( rst )=valid_measures

require(survival)
rst = NULL;
for (m in S4_valid_measures){
	metabolite = S4[, m]
	model = coxph(Surv(apo_time, inz_apo) ~  log(metabolite) +
					ltalteru + as.factor(lcsex) + ltbmi## model 1
			#+(lp_diab_who06==4|lp_diab_who06==5)  ##model 2
			#+log(ll_choln)+log(ll_hdln)+log(ltsysmm)+ as.factor(ltcigreg)##model 3
			#+total2HDL ##model 4
			#+ltdiamm ## model 5
			,subset = which(S4$prev_apo == 0&!(S4$apo_typ %in% c(1,2,3,4,9))),
			S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = "Stroke survival analysis_model1.csv")

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


###################	Myocardial vascualr disease ##############################
require(survival)
rst = NULL;
for (m in S4_valid_measures){
	metabolite = S4[, m]
	model = coxph(Surv(apo_time, inz_mi) ~  log(metabolite) +
					ltalteru + as.factor(lcsex) + ltbmi## model 1
					#+(lp_diab_who06==4|lp_diab_who06==5)  ##model 2
					#+log(ll_choln)+log(ll_hdln)+log(ltsysmm)+ as.factor(ltcigreg)##model 3
	 				#+total2HDL ##model 4
					#+ltdiamm ## model 5
					,subset = which(S4$prev_mi == 0),
					S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = "MI survival analysis_model1.csv")

plot(survfit(Surv(mi_time, S4$inz_mi)~(log(S4$PC_aa_C32_2) > 1.2), S4, subset= which(S4$prev_mi == 0)), log = "y", col = c("red","green"))

###variable selection by boosting method
metabo.asso = scan(what = character())
C14_1
C16_1
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
PC_aa_C34_4
PC_aa_C36_2
PC_aa_C36_3
PC_aa_C36_6
PC_ae_C36_2
PC_ae_C38_0
PC_ae_C38_2
PC_ae_C40_1
SM_C24_1


clinical = c("ltalteru", "ltbmi", "lcsex","lp_diab_who06", "ltsysmm", "ll_hdln", "ll_choln", "ltcigreg")

#### stepwise selection of cox regression
selectCox <- function(formula, data, rule = "aic") {
	require("rms")
	require("prodlim")
	fit <- cph(formula, data, surv = TRUE)
	bwfit <- fastbw(fit, rule = rule)
	if (length(bwfit$names.kept) == 0) {
		newform <- reformulate("1", formula[[2]])
		newfit <- prodlim(newform, data = data)
	} else{
		newform <- reformulate(bwfit$names.kept, formula[[2]])
		newfit <- cph(newform, data, surv = TRUE)
	}
	out <- list(fit = newfit,In = bwfit$names.kept)
	out$call <- match.call()
	class(out) <- "selectCox"
	out
}

selectCox(
		formula = Surv(mi_time, inz_mi) ~  Arg+Trp+lysoPC_a_C16_0+lysoPC_a_C16_1+lysoPC_a_C17_0+lysoPC_a_C18_0+lysoPC_a_C18_1+lysoPC_a_C18_2+PC_aa_C28_1+PC_aa_C32_2+PC_aa_C34_2+PC_aa_C34_3+PC_aa_C34_4+PC_aa_C36_2+PC_aa_C36_3+PC_aa_C36_6+PC_ae_C36_2+PC_ae_C38_0+PC_ae_C40_1+SM_C24_1,
		data = data.frame(log(S4[, metabolites.asso]), S4[, model1], mi_time = S4$mi_time, inz_mi = S4$inz_mi)[subset, ],
		rule = "p")

+ ltalteru + log(ltdiamm) + log(ltsysmm) + log(ll_hdln) + log(ll_choln) + as.factor(lp_diab_who06) + as.factor(lcsex) + as.factor(ltcigreg)


