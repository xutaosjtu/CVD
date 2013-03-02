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

#########################	hypertension	################################
#Prediction of future hypertension
table(S4$lthyact, S4$uthyact, useNA='ifany')
rst = NULL;
for (m in S4_valid_measures){
	metabolite = S4[, m]
	model = glm(as.factor(uthyact) ~  metabolite +
					ltalteru + as.factor(lcsex) + ltbmi## model 1
					+ ltsysmm #+ ltdiamm ##model 2
			+(lp_diab_who06==4|lp_diab_who06==5)  ##model 2
			+log(ll_chola)+log(ll_hdla)+ total2HDL##model 3
			+ my.cigreg + my.alkkon ##model 4
			+as.factor(ltantihy) #+ ltmata #model 6 medication
			,subset = which(S4$lthyact == 2),
			S4, family = binomial(link = "logit"))
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = "Hypertension survival analysis_model2.csv")


#Association with current hypertension
S4.feature = scan(what = character())
ltalteru
lcsex
ltbmi
my.diab
ll_chola
ll_hdla
ll_ldla
total2HDL
my.cigreg
ltalkkon
ltsysmm
ltdiamm
ltmbbl
ltmace
ltmata
ltantihy
my.alkkon

F4.feature = scan(what = character())
utalteru
ucsex
utbmi
my.diab
ul_chola
ul_hdla
ul_ldla
total2HDL
my.cigreg
utalkkon
utsysmm
utdiamm
utmbbl
utmace
utmata
utantihy
my.alkkon

#F4 = preprocess(F4, F4_valid_measures)
require(nlme)
valid_measures = intersect(S4_valid_measures, F4_valid_measures)
participants=rep(1:1009,2)
tmpF4 = F4[Cohort$zz_nr_f4,F4.feature]
colnames(tmpF4) = S4.feature
data=data.frame(
		participants, 
		disease =  as.factor(2-c(S4[Cohort$zz_nr_s4, "lthyact"], F4[Cohort$zz_nr_f4, "uthyact"])),
		rbind(S4[Cohort$zz_nr_s4, S4.feature], tmpF4)
)
data$platform = rep(1:2, each = 1009)
data$ltantihy = 2-data$ltantihy

sub = which(S4[Cohort$zz_nr_s4,"ltantihy"]==F4[Cohort$zz_nr_f4,"utantihy"])
sub = which(S4[Cohort$zz_nr_s4,"lthyact"]==1&F4[Cohort$zz_nr_f4,"uthyact"]==1)
sub = which(S4[Cohort$zz_nr_s4,"lthyact"]!=F4[Cohort$zz_nr_f4,"uthyact"])
sub = c(sub, sub+1009)

## medication in the populetaion at baseline and follow-up
table(S4[Cohort$zz_nr_s4[subset],"ltantihy"],F4[Cohort$zz_nr_f4[subset],"utantihy"], S4[Cohort$zz_nr_s4[subset],"lthyact"], F4[Cohort$zz_nr_f4[subset],"uthyact"])


rst=NULL
for(i in valid_measures){
	data$m = c(log(S4[Cohort$zz_nr_s4, i]), log(F4[Cohort$zz_nr_f4, i]))
	#data$m = scale(data$m)
#	mixed.dum <- lme( m ~ disease + 
#					as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
#					platform +
#					ltalteru + as.factor(lcsex) + ltbmi## model 1
#					+ my.cigreg + my.alkkon ##model 2
#					+ my.diab  ##model 3
#					+ ll_chola+ll_hdla##model 4
#			,random = ~  1 | participants, na.action=na.exclude, 
#			#subset = sub,
#			data=data[sub,])
#	rst = rbind(rst, summary(mixed.dum)$tTable[2,])
	mixed.dum <- glmer(disease ~ m + 
					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
					platform +
					ltalteru + as.factor(lcsex) + ltbmi## model 1
					+ my.cigreg + my.alkkon ##model 2
					+ my.diab  ##model 3
					+ ll_chola+ll_hdla##model 4
			        + (1 | participants), 
			na.action=na.exclude, family = binomial,
			data)
	rst = rbind(rst, summary(mixed.dum)@coefs[2,])
}
rst=data.frame(rst,
		fdr=p.adjust(rst[, 4], method="fdr"),
		bonf=p.adjust(rst[, 4], method="bonferroni")
)
rownames(rst)=valid_measures
write.csv(rst, file = "hypertension associated metabolites_full model_logistic.csv")

plot(res~fit,data=data.frame(fit = mixed.dum$fitted[,1],res = mixed.dum$residuals[,1]))
abline(lm(res~fit,data=data.frame(fit = mixed.dum$fitted[,1],res = mixed.dum$residuals[,1])))

############	calculate the residues	########
data = data.frame(
		"alteru" = S4$ltalteru, 
		"sex" = as.factor(S4$lcsex),
		"bmi" = S4$ltbmi,
		"cigreg" = as.factor(S4$ltcigreg),
		"alkkon" = as.factor(S4$my.alkkon),
		"diabetes" = as.factor((S4$lp_diab_who06==4|S4$lp_diab_who06==5)),
		"chol" = log(S4$ll_chola),
		"HDL" = log(S4$ll_hdla),
		#"med1" = as.factor(S4$ltmbbl),
		#"med2" = as.factor(S4$ltmace),
		#"med3" = as.factor(S4$ltmata),
		'med' = as.factor(S4$ltantihy),
		"hyper" = S4$lthyact,
		log(S4[,valid_measures])
)
data = data[Cohort$zz_nr_s4, ]
S4.residue = residue(data, valid_measures, adj = colnames(data)[1:12], control_group = 1:dim(data)[1])

data = data.frame(
		"alteru" = F4$utalteru, 
		"sex" = as.factor(F4$ucsex),
		"bmi" = F4$utbmi,
		"cigreg" = as.factor(F4$utcigreg),
		"alkkon" = as.factor(F4$my.alkkon),
		"diabetes" = as.factor((F4$uk_diab_who06==4|F4$uk_diab_who06==5)),
		"chol" = log(F4$ul_chola),
		"HDL" = log(F4$ul_hdla),
		#"med1" = as.factor(F4$utmbbl),
		#"med2" = as.factor(F4$utmace),
		#"med3" = as.factor(F4$utmata),
		'med' = as.factor(F4$utantihy),
		"hyper" = F4$uthyact,
		log(F4[,valid_measures])
)
data = data[Cohort$zz_nr_f4, ]
F4.residue = residue(data, valid_measures, adj = colnames(data)[1:10], control_group = 1:dim(data)[1])

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

############################	Stroke	########################################
#require(survival)
#rst = NULL;
#for (m in S4_valid_measures){
#	metabolite = S4[, m]
#	model = coxph(Surv(apo_time, inz_apo) ~  log(metabolite) +
#					ltalteru + as.factor(lcsex) + ltbmi## model 1
#			+(lp_diab_who06==4|lp_diab_who06==5)  ##model 2
#			+log(ll_choln)+log(ll_hdln)+log(ltsysmm)+ as.factor(ltcigreg) + ltalkkon ##model 3
#	 		#+ log(lh_crp) + total2HDL ##model 4
#			#+ltdiamm ## model 5
#			,subset = which(S4$prev_apo == 0&!(S4$apo_typ %in% c(5,1,2,9))),#
#			S4)
#	rst = rbind(rst, summary(model)$coefficients[1,])
#}
#rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
#rownames(rst) = S4_valid_measures
#write.csv(rst, file = "Stroke survival analysis_model3.csv")
#
#table(S4$inz_apo==1, S4$apo_typ)
#S4$my.apo_typ = S4$apo_typ
#S4$my.apo_typ[which(S4$apo_typ == 0)] = 0
#S4$my.apo_typ[which(S4$apo_typ == 2)] = 1
#S4$my.apo_typ[which(S4$apo_typ == 3)] = 2
#S4$my.apo_typ[which(S4$apo_typ == 4)] = 2
#S4$my.apo_typ[which(S4$apo_typ == 5)] = 3
#S4$my.apo_typ[which(S4$apo_typ == 9)] = 4
#
#require(caret)
#require(pls)
#require(gplots)
#index = which((S4$my.apo_typ == 1|S4$my.apo_typ == 2)&S4$prev_apo==0)#
#S4pls<-plsr(S4$my.apo_typ[index] ~ . , data=log2(S4[index, c(S4_valid_measures)]),  validation = "CV")
#S4pls<-plsda(x=log2(S4[index, c(S4_valid_measures)]), y = as.factor(S4$my.apo_typ[index]))
##S4pls<-plsda(x=data.normalized, COPD.data$COPD, ncomp = 10)
#
#color = greenred(24)[c(c(4:1)*2,c(18,20, 22, 24)) ]
#plot(S4pls$scores[,c(1,2)],col = color[c(1,8)][S4$my.apo_typ[index]], pch=c(17, 19)[S4$my.apo_typ[index]])
#legend(1, -3, 
#		legend = c("Ischemic","Hemorrhagic"), 
#		col = c(1:8), pch = c(17,19)
#)
#
#S4pca = prcomp(log(S4[index, S4_valid_measures]) )
#S4pca = pcr(S4$my.apo_typ[index] ~ ., data = log2(S4[index, c(S4_valid_measures)]))
#plot(S4pca$scores,col = color[S4$my.apo_typ[which(S4$my.apo_typ == 1|S4$my.apo_typ == 2)]], pch=c(17, 19)[S4$my.apo_typ[which(S4$my.apo_typ == 1|S4$my.apo_typ == 2)]])
#
#difference = scan(what = character())
#PC_aa_C28_1
#PC_ae_C38_1
#SM__OH__C22_1
#SM_C16_0
#SM_C24_0
#C0
#C10_2
#Pro
#Taurine
#PC_ae_C40_4


###################	Myocardial infarction ##############################
S4$my.cigreg = S4$ltcigreg
S4$my.cigreg[which(S4$ltcigreg ==2)] = 1
S4$my.cigreg[which(S4$ltcigreg ==3)] = 2
S4$my.cigreg[which(S4$ltcigreg ==4)] = 3
S4$my.cigreg = 3-S4$my.cigreg
S4$my.cigreg = factor(S4$my.cigreg, ordered = F)

S4$total2HDL = S4$ll_chola/S4$ll_hdla

require(survival)
rst = NULL; rst1 = NULL
rst2 = NULL; rst3 = NULL
for (m in metabo.ratio.asso){
	metabolite = scale((S4[, m]))
	model = coxph(Surv(mi_time, inz_mi) ~ metabolite +
					scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
					#+ my.diab  ##model 2
					#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
	 				#+ scale(lh_crp)  ##model 4
					#+ as.factor(ltmstati)& S4$ltmstati !=1
					,subset = which(S4$prev_mi == 0 ),
					S4)
	rst = rbind(rst, summary(model)$coefficients[1,])
	#rst1 = rbind(rst , summary(model)$coefficients[10,])
	#rst2 = rbind(rst2, summary(model)$coefficients[11,])
	#rst3 = rbind(rst3, summary(model)$coefficients[12,])
	#table(model$y[,2]) #number of sample used exactly in the estimation.
}
table(model$y[,2])
rst = data.frame(rst, FDR = p.adjust(rst[,5], method = "BH"), bonferroni = p.adjust(rst[,5], method = "bonferroni"))
rownames(rst) = metabo.ratio.asso
#rst = cbind(rst, annotation[rownames(rst),])
write.csv(rst, file = "metabolites ratio invest_MI survival analysis_model1.csv")

plot(survfit(Surv(mi_time, S4$inz_mi)~(log(S4$PC_aa_C32_2) > 1.2), S4, subset= which(S4$prev_mi == 0)), log = "y", col = c("red","green"))

############	Hazardous ratios in different quantiles	################
require(gplots)
pdf("quintile plot of replative risk (ynorm) _full model _decile.pdf", width = 12, height = 12)
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
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2
PC_aa_C32_3
PC_aa_C34_2
PC_aa_C34_3
PC_aa_C36_2
PC_aa_C36_3
PC_ae_C36_1
PC_ae_C36_2
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



