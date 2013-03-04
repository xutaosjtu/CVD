# regression analysis
#	model 1: basic model
	feature.cont = c("ltalteru", "ltalkkon" , "ltbmi")
 	feature.disc = c("lcsex")
#	model 2: full model
	feature.cont = 
	feature.disc = c("lcsex")
#	model 3: Estimated clinical model
#	Hypertension:
	feature.cont = c("ltalteru","ll_hgb", "ll_trin", "ltbmi" , "lcsex" , "lh_crp", "ll_hdln") 
	feature.disc = c("lcsex")
#	Stroke: 
	feature.cont = c("ltalteru" , "lh_crp" , "ll_choln" )
	feature.disc = c("lcsex")
#	Myocardial infarction: 
	feature.cont = c("ll_choln" , "ll_ldln" , "lthumf" )
	feature.disc = c("lcsex")
#	Heart failure
	feature.cont = c("ltalteru"	, "ll_choln" , "ll_hdln" , "ll_trin" , "lthumf" , "lttumf" , "ltalkkon" , "ll_hgb")
	feature.disc = c("lcsex")
#	model 4: CVD scoring system
#	Framingham Heart Study:
	feature.cont = c("ltalteru", "ll_hdln" , "ll_choln" , "ltsysmm")
	feature.disc = c("lcsex" , "lp_diab_who06")
#	The SCORE project:
	feature.cont = c("ltalteru", "ltsysmm" , "ltotal2HDL")
	feature.disc = c("lcsex")
# 	Hypertension:
	feature.cont = c("ltbmi", "ltalteru", "ltalkkon", "lttumf", "waist2hip")
	feature.disc = c("ltrauchp", "lcsex")
# Author: tao.xu
###############################################################################



# Estimation of clinical model:
data = data.frame (disease = S4[,diseases[4]], S4[,feature.cont], as.factor(S4[,feature.disc]))
data <- na.omit(data)
model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease != 3) )
smodel = step(model, trace = 0)

data = data.frame (disease = F4[,diseases[4]], F4[,feature.cont], as.factor(F4[,feature.disc]))
data <- na.omit(data)
model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease != 3))
smodel = step(model, trace = 0)

#, subset = c(which(data$disease==1),which(data$disease==4))
######################	S4	###############

rst = NULL
for (m in S4_valid_measures){
	data = data.frame (disease = as.factor(S4[,diseases.s4[1]]), metabolite = S4[,m], S4[,feature.cont], apply(S4[,feature.disc], 2, function(x) as.factor(x)))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2:4,4] )
}
rst = data.frame(rst, dis1_FDR = p.adjust(rst[,1], method = "BH"), dis2_FDR = p.adjust(rst[,2], method = "BH"), dis3_FDR = p.adjust(rst[,3], method = "BH"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = paste( names(diseases.s4)[1], "categorical model_4 S4.csv", sep = " "))

rst = NULL
for (m in S4_valid_measures){
	data = data.frame (disease = S4[,diseases.s4[1]], metabolite = log(S4[,m]), S4[,feature.cont], apply(S4[,feature.disc], 2, function(x) as.factor(x)))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2,] )
}
rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = paste( names(diseases.s4)[1], "model_3 S4.csv", sep = " "))

for(i in 4:4){
	rst = NULL
	for(m in S4_valid_measures){
		data = data.frame (disease = S4[,diseases[i]], metabolite = log(S4[,m]), S4[,feature.cont], as.factor(S4[,feature.disc]))
		
		model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease!=3))
		rst = rbind ( rst , summary(model)$coefficients[2,] )
	}
	rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
	rownames(rst) = S4_valid_measures
	write.csv(rst, file = paste( names(diseases)[i], "model_3 S4.csv", sep = " " ))
}	

rst = NULL
for(m in S4_valid_measures){
	data = data.frame (disease = S4[,diseases[3]], S4[,m], S4[,feature.cont], apply(S4[,feature.disc], 2, function(x) as.factor(x)))
	
	model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease!=3), na.action = na.exclude)
	
	rst = rbind ( rst , summary(model)$coefficients[2,] )
}


######################	F4	###############
# Hypertension
rst = NULL
for (m in F4_valid_measures){
	data = data.frame (disease = as.factor(F4[,diseases.f4[1]]), metabolite = F4[,m], F4[,feature.cont], apply(F4[,feature.disc], 2, function(x) as.factor(x)))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2:4,4] )
}
rst = data.frame(rst, dis1_FDR = p.adjust(rst[,1], method = "BH"), dis2_FDR = p.adjust(rst[,2], method = "BH"), dis3_FDR = p.adjust(rst[,3], method = "BH"))
rownames(rst) = F4_valid_measures
write.csv(rst, file = paste( names(diseases.f4)[1], "categorical model_4 F4.csv", sep = " "))

rst = NULL
for (m in F4_valid_measures){
	data = data.frame (disease = F4[,diseases.f4[1]], metabolite = log(F4[,m]), F4[,feature.cont], apply(F4[,feature.disc], 2, function(x) as.factor(x)))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2,] )
}
rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
rownames(rst) = F4_valid_measures
write.csv(rst, file = paste( names(diseases.f4)[1], "model_4 F4.csv", sep = " "))

# stroke, myocardio infarction, heart failure
for(i in 4:4){
	rst = NULL
	for(m in F4_valid_measures){
		data = data.frame (disease = F4[,diseases.f4[i]], metabolite = log(F4[,m]), F4[,feature.cont], as.factor(F4[,feature.disc]))
		
		model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease!=3))
		rst = rbind ( rst , summary(model)$coefficients[2,] )
	}
	rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
	rownames(rst) = F4_valid_measures
	write.csv(rst, file = paste( names(diseases.f4)[i], "model_3 F4.csv", sep = " " ))
}	

################# Combining three cardiovascular diseases () together
tmp = (S4[,diseases.s4[2:4]])
D = apply(tmp, 1, function(x) sum(x==1))
table(D)
S4$CVD =  apply(tmp, 1, function(x) sum(x==1, na.rm = T))!=0

tmp = (F4[,diseases.s4[2:4]])
table(apply(tmp, 1, function(x) sum(x==1)))
F4$CVD = apply(tmp, 1, function(x) sum(x ==1, na.rm = T)!=0)

tmp = (S4[,diseases.s4[2:4]])
table(apply(tmp, 1, function(x) sum(x==1)))

#	regression analysis of S4]
substr(feature.cont, 1, 2) <- "l"
substr(feature.disc, 1, 2) <- "l"
rst = logisticRegression(S4, S4$CVD, S4_valid_measures, feature.cont, feature.disc)
write.csv(rst, file = "CVD combined_Model Framingham_S4.csv")

# regression analysis of F4
substr(feature.cont, 1, 2) <- "u"
substr(feature.disc, 1, 2) <- "u"
rst = logisticRegression(F4, F4$CVD, F4_valid_measures, feature.cont, feature.disc)
write.csv(rst, file = "CVD combined_Model SCORE_F4.csv")

#################	Ratio analysis
concen2ratio = function(data){
	
	rst = NULL
	
	for(i in 1:dim(data)[2]){
		tmp =  t(apply(data, 1, function(x) x[i]/x[-i]))
		colnames(tmp) = paste(colnames(data)[i], colnames(data)[-i], sep = ".")
		rst = cbind(rst, tmp)
	
	}
	
	return(rst)

}

substr(feature.cont, 1, 2) <- "l"
substr(feature.disc, 1, 2) <- "l"
S4_ratio = concen2ratio(S4[, S4_valid_measures[c(1:21)]])
rst = logisticRegression(data.frame(S4 , S4_ratio), S4$lc044f_1, colnames(S4_ratio), feature.cont, feature.disc, metalog = FALSE)
colnames(S4_ratio)[which(rst[,5]<0.05)]


substr(feature.cont, 1, 2) <- "u"
substr(feature.disc, 1, 2) <- "u"
F4_ratio = concen2ratio(F4[, F4_valid_measures[c(1:23)]])
rst = logisticRegression(data.frame(F4 , F4_ratio) , F4$us_c04a, colnames(F4_ratio), feature.cont, feature.disc, metalog = FALSE)
colnames(F4_ratio)[which(rst[,5]<0.05)]


## test the association of the three biomarkers in the cross-sectional data
## F4
rst = NULL;
for(m in F4_valid_measures){
	F4$m = log(F4[,m])
	model = glm(as.factor(utmi)~ m 
					+ scale(utalteru) + as.factor(ucsex) + scale(utbmi)## model 1
					+ my.diab  ##model 2
					+ scale(utsysmm) + my.cigreg + my.alkkon  + scale(ul_chola) + scale(ul_hdla) ##model 3+ total2HDL
					+ scale(uh_crp)  ##model 4
					#+ as.factor(my.medication),
			,family=binomial,
			F4)
#	model = clogit((utmi-1) ~ m 
#					+ scale(utalteru) + as.factor(ucsex) + scale(utbmi)## model 1
#					+ my.diab  ##model 2
#					+ scale(utsysmm) + my.cigreg + my.alkkon  + scale(ul_chola) + scale(ul_hdla) ##model 3+ total2HDL
#					+ scale(uh_crp)  ##model 4
#					+strata(my.medication)
#			,data = F4)
	rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst)=F4_valid_measures


## S4
rst = NULL;
for(m in S4_valid_measures){
	S4$m = log(S4[,m])
	model = glm(as.factor(ltmi)~m 
					+ scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
					+ my.diab  ##model 2
					+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_ldla) ##model 3+ total2HDL
					+ scale(lh_crp)  ##model 4
					#+my.medication
	,family =binomial
	#,subset =which(S4$my.medication)
	, S4)
#	model = clogit(prev_mi~m 
#					+ scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
#					+ my.diab  ##model 2
#					+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_ldla) ##model 3+ total2HDL
#					+ scale(lh_crp)  ##model 4
#					+strata(ltmi)
#	,data = S4)
	rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst)=S4_valid_measures

s4.rst = read.csv("MI associated metabolites_S4_logistic.csv")
f4.rst = read.csv("MI associated metabolites_F4_logistic.csv")
combined.rst=merge(s4.rst, f4.rst, by="X", sort=F, all=T)
write.csv(combined.rst, file = "MI associated metabolties_cross sectional_logistic.csv", quote=F, row.names=F)



C0
C4
C5
C4_1_DC___C6
C8
C8_1
C10
C10_1
C10_2
C14_2
C16
C16_1
C18
C18_1
C18_2
Phe
PC_aa_C36_4
PC_aa_C38_4
PC_ae_C32_1
PC_ae_C34_0
PC_ae_C34_3
PC_ae_C36_2
PC_ae_C38_2
PC_ae_C40_2
PC_ae_C40_5
PC_ae_C42_0
PC_ae_C42_2
PC_ae_C42_4
PC_ae_C42_5
lysoPC_a_C20_4
SM__OH__C14_1









