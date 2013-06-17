# The main focuse of the analysis rest in two parts:
# 1. To find better predictors other than the conventional biomarkers of hypertension
# 2. To find metabolites associated with existing hypertension
# Author: tao.xu
###############################################################################


S4.feature = scan(what = character())
ltalteru
lcsex
ltbmi
my.diab
ll_chola
ll_hdla
ll_ldla
lh_crp
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
my.hyper
lthyact

F4.feature = scan(what = character())
utalteru
ucsex
utbmi
my.diab
ul_chola
ul_hdla
ul_ldla
uh_crp
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
my.hyper
uthyact

print("medication in the population at baseline and follow-up")
table(S4[Cohort$zz_nr_s4,"ltantihy"],F4[Cohort$zz_nr_f4,"utantihy"], S4[Cohort$zz_nr_s4,"lthyact"], F4[Cohort$zz_nr_f4,"uthyact"])

## Prediction of future hypertension
print("hypertension status at baseline and follow-up")
table('S4' = S4$lthyact, 'F4' = S4$uthyact, useNA='ifany')
rst = NULL;
for (m in valid_measures){
	## fixed effect estimation
	metabolite = scale(log(S4[, m]))
	model = glm(as.factor(2-uthyact) ~  metabolite +
					ltalteru + as.factor(lcsex) # crude model
					+ scale(ltsysmm) # multivariate model
					+ scale(ltbmi)
					+ as.factor(my.cigreg) + as.factor(my.alkkon)
					+ my.diab
					+ scale(ll_chola) + scale(ll_hdla)
					+ scale(log(lh_crp))
			,subset = which(S4$lthyact == 2),
			S4, family = binomial(link = "logit"))
	rst = rbind(rst, summary(model)$coefficients[5,])
	
	## mixed effect estimation
#	data$m = c(log(S4[Cohort$zz_nr_s4, m]), log(F4[Cohort$zz_nr_f4, m]))
#	mixed.dum <- glmer(disease ~ m + 
#					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
#					as.factor(platform) +
#					ltalteru + as.factor(lcsex) ## model 1
#					+ ltsysmm
#					+ ltbmi + my.cigreg + my.alkkon ##model 2
#					+ my.diab  ##model 3
#					+ ll_chola+ll_hdla##model 4
#					+ (1 | participants), 
#			na.action=na.exclude, family = binomial,
#			data[sub,])
#	rst = rbind(rst, summary(mixed.dum)@coefs[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = valid_measures
write.csv(rst, file = "Hypertension survival analysis_mixed effect_full model.csv")


## Association with current hypertension
# fixed effect model: logistic regression
rst = NULL # association in S4 
for (m in valid_measures){
	metabolite = scale(log(S4[, m]))
	model = glm(as.factor(2-lthyact) ~  metabolite +
					ltalteru + as.factor(lcsex) # basic model
#					+ scale(ltbmi) # multivariate model
#					+ as.factor(my.cigreg) + as.factor(my.alkkon)
#					+ my.diab
#					+ scale(ll_chola) + scale(ll_hdla)
#					+ scale(log(lh_crp))
			,#subset = which(S4$ltantihy == 2), # antihypertensive medication 
			S4, family = binomial(link = "logit"))
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = valid_measures
write.csv(rst, file = "Hypertension cross-sectional_S4_include hyper(med)_crude model_normed.csv")

rst = NULL # association in F4
for (m in valid_measures){
	metabolite = scale(log(F4[, m]))
	model = glm(as.factor(2-uthyact) ~  metabolite +
					utalteru + as.factor(ucsex) ## model 1
					+ scale(utbmi) # multivariate model
					+ as.factor(my.cigreg) + as.factor(my.alkkon)
					+ my.diab
					+ scale(ul_chola) + scale(ul_hdla)
					+ scale(log(uh_crp))
			,#subset = which(F4$utantihy == 2), # antihypertensive medication
			F4, family = binomial(link = "logit"))
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = valid_measures
write.csv(rst, file = "Hypertension cross-sectional_F4_include hyper(med)_full model.csv")


valid_measures = intersect(S4_valid_measures, F4_valid_measures)
participants=rep(1:length(Cohort$zz_nr_s4),2)
tmpF4 = F4[Cohort$zz_nr_f4,F4.feature]
colnames(tmpF4) = S4.feature
data=data.frame(
		participants, 
		disease =  as.factor(2-c(S4[Cohort$zz_nr_s4, "lthyact"], F4[Cohort$zz_nr_f4, "uthyact"])),
		rbind(S4[Cohort$zz_nr_s4, S4.feature], tmpF4)
)
data$platform = rep(1:2, each = length(Cohort$zz_nr_s4))
data$ltantihy = 2-data$ltantihy

#sub = which(S4[Cohort$zz_nr_s4,"ltantihy"]==F4[Cohort$zz_nr_f4,"utantihy"])
#sub = which(S4[Cohort$zz_nr_s4,"lthyact"]==1&F4[Cohort$zz_nr_f4,"uthyact"]==1)
#sub = which(S4[Cohort$zz_nr_s4,"lthyact"]!=F4[Cohort$zz_nr_f4,"uthyact"])
#sub = which(S4[Cohort$zz_nr_s4,"lthyact"]==2)
#
#sub = c(sub, sub+1009)

# mixed effect model: linear mixed model
require(nlme)
rst=NULL
for(i in valid_measures){
	data$m = c(scale(log(S4[Cohort$zz_nr_s4, i])), scale(log(F4[Cohort$zz_nr_f4, i])))
	#data$m = scale(data$m)
	mixed.dum <- lme( m ~ disease + #linear mixed model
					#ltdiamm + ltantihy +
					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model medication	
					#platform +
					ltalteru + as.factor(lcsex) # crude model
					+ ltbmi+ my.cigreg + my.alkkon + my.diab + ll_chola+ll_hdla +log(lh_crp) # multivariate model
			,random = ~  1 | participants, na.action=na.exclude, 
			subset = which(data$ltantihy!=1),
			data=data[order(data$participants),])
	rst = rbind(rst, summary(mixed.dum)$tTable[2,])
#	mixed.dum <- glmer(disease ~ m + # logisitic mixed model
#					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
#					platform +
#					ltalteru + as.factor(lcsex) + ltbmi # crude model
#					+ my.cigreg + my.alkkon # multivariate model
#					+ my.diab
#					+ ll_chola+ll_hdla
#			        + (1 | participants), 
#			na.action=na.exclude, family = binomial,
#			data)
#	rst = rbind(rst, summary(mixed.dum)@coefs[2,])
}
rst=data.frame(rst,
		fdr=p.adjust(rst[, 5], method="fdr"),
		bonf=p.adjust(rst[, 5], method="bonferroni")
)
rownames(rst)=valid_measures
write.csv(rst, file = "Hypertension associated metabolites_without medication_full model.csv")

plot(res~fit,data=data.frame(fit = mixed.dum$fitted[,1],res = mixed.dum$residuals[,1]))
abline(lm(res~fit,data=data.frame(fit = mixed.dum$fitted[,1],res = mixed.dum$residuals[,1])))

# General estimate equation
require(gee) 
#require(geepack)
rst=NULL 
for(i in valid_measures){
	data$m = c(scale(log(S4[Cohort$zz_nr_s4, i])), scale(log(F4[Cohort$zz_nr_f4, i])))#normalize to comparable value
	tmp = data[order(data$participants), ]
	model <- gee(disease ~ m +
					#ltdiamm + ltantihy +
					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
					#platform +
					ltalteru + as.factor(lcsex) ## model 1
					+ ltbmi+ my.cigreg + my.alkkon + my.diab + ll_chola+ll_hdla +log(lh_crp) ##model 4
			,id = participants,# na.action=na.exclude, 
			data= tmp,
			subset = which(!is.na(tmp$disease)),#
			corstr = "exchangeable",
			family = binomial)
	rst = rbind(rst, summary(model)$coef[2,])
}
rst = data.frame(rst, pvalue = 2*pnorm(-abs(rst[,5])))
rst=data.frame(rst,
		fdr=p.adjust(rst$pvalue, method="fdr"),
		bonf=p.adjust(rst$pvalue, method="bonferroni")
)
rownames(rst)=valid_measures
write.csv(rst, file = "Hypertension associated metabolites_crude model_GEE (without med).csv")

## Association of blood pressure with metabolite concentrations
rst=NULL; #linear mixed model
for(i in valid_measures){
	data$m = c(scale(S4[Cohort$zz_nr_s4, i]), (scale(F4[Cohort$zz_nr_f4, i])))
	mixed.dum <- lme( ltdiamm ~ m + #linear mixed model
					#ltdiamm + ltantihy +
					#as.factor(ltantihy) + #+as.factor(ltmbbl) #model medication	
					#platform +
					ltalteru + as.factor(lcsex) # crude model
					+ ltbmi+ my.cigreg + my.alkkon + my.diab + ll_chola+ll_hdla +log(lh_crp) # multivariate model
			,random = ~  1 | participants, na.action=na.exclude, # multivariate model 
			subset = which(data$ltantihy!=1),
			data=data)
	rst = rbind(rst, summary(mixed.dum)$tTable[2,])
}
rst=data.frame(rst,
		fdr=p.adjust(rst[, 5], method="fdr"),
		bonf=p.adjust(rst[, 5], method="bonferroni")
)
rownames(rst)=valid_measures
write.csv(rst, file = "DiastolicBP associated metabolites_without medication_full model_unnorm_lme.csv")

rst=NULL; #GEE model
for(i in valid_measures){
  data$m = c(scale(log(S4[Cohort$zz_nr_s4, i])), scale(log(F4[Cohort$zz_nr_f4, i])))#normalize to comparable value
  tmp = data[order(data$participants), ]
  model <- gee(ltsysmm ~ m +
                 #ltdiamm + ltantihy +
                 #as.factor(ltantihy) + #+as.factor(ltmbbl) #model 5 medication	
                 #platform +
                 ltalteru + as.factor(lcsex) ## model 1
               + ltbmi+ my.cigreg + my.alkkon + my.diab + ll_chola+ll_hdla +log(lh_crp) ##model 4
               ,id = participants,# na.action=na.exclude, 
               data= tmp,
               subset = which(!is.na(tmp$disease)&tmp$ltantihy!=1),#
               corstr = "exchangeable")
  rst = rbind(rst, summary(model)$coef[2,])
}
rst = data.frame(rst, pvalue = 2*pnorm(-abs(rst[,5])))
rst=data.frame(rst,
               fdr=p.adjust(rst$pvalue, method="fdr"),
               bonf=p.adjust(rst$pvalue, method="bonferroni")
)
rownames(rst)=valid_measures
write.csv(rst, file = "DiastolicBP associated metabolites_without medication_crude model_norm_GEE.csv")

## Association with current hypertension
# cross-sectional linear regression
rst = NULL # association in S4 
for (m in valid_measures){
  metabolite = scale(log(S4[, m]))
  model = lm(ltsysmm ~  metabolite +
                ltalteru + as.factor(lcsex) # basic model
      					+ scale(ltbmi) # multivariate model
      					+ as.factor(my.cigreg) + as.factor(my.alkkon)
      					+ my.diab
      					+ scale(ll_chola) + scale(ll_hdla)
      					+ scale(log(lh_crp))
             #+as.factor(ltantihy)
              ,subset = which(S4$ltantihy == 2), # antihypertensive medication 
              S4)
  rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = valid_measures
write.csv(rst, file = "distolicBP cross-sectional_S4_include hyper(med)_crude model_normed.csv")

rst = NULL # association in F4
for (m in valid_measures){
  metabolite = scale(log(F4[, m]))
  model = lm(utsysmm ~  metabolite +
                utalteru + as.factor(ucsex) ## model 1
              #+ scale(utbmi) # multivariate model
              #+ as.factor(my.cigreg) + as.factor(my.alkkon)
              #+ my.diab
              #+ scale(ul_chola) + scale(ul_hdla)
             # + scale(log(uh_crp))
              ,#subset = which(F4$utantihy == 2), # antihypertensive medication
              F4, family = binomial(link = "logit"))
  rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = valid_measures
write.csv(rst, file = "Hypertension cross-sectional_F4_include hyper(med)_full model.csv")

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
#close.screen(all = T)
par(mfrow = c(6,3))
k = 1
for(i in metabo.selected){
	#screen(k)
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


