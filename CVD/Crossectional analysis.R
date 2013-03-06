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

## By a matched case control study
data = S4[which(S4$prev_mi==0),c("ltalteru", "lcsex", "ltbmi","my.diab","ltsysmm","my.cigreg","my.alkkon","ll_chola", "ll_hdla", "ll_ldla", "lh_crp", "inz_mi", S4_valid_measures)]
data = na.omit(data)
data$my.cigreg = as.numeric(data$my.cigreg)
data$my.diab = as.numeric(data$my.diab)
data$Arg_Trp = data$Arg/data$Trp

#Y = 
tr = data$inz_mi
x = data[,c("ltalteru", "lcsex", "my.diab","my.cigreg","my.alkkon")]
rst.match = Match(Tr = tr, X =x, exact = T, M = 1, ties = F)

#rst.balance = MatchBalance(inz_mi ~ scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
#+ my.diab  ##model 2
#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_ldla) ##model 3+ total2HDL
#+ scale(lh_crp)  ##model 4
#, data, match = rst.match, nboots = 10)

data = data[c(rst.match$index.control, unique(rst.match$index.treated)), ]
#data$ID = rep(1:table(data$inz_mi)[1], 2)
data$ID = c(rep(1:table(data$inz_mi)[2], each = 1),1:table(data$inz_mi)[2])
require(survival)
rst = NULL
for(m in S4_valid_measures){
	#data$m = log(data[,"Arg"])-log(data$Trp)
	data$m = log(data[,m])
	model = clogit(inz_mi ~ m 
	          		+ scale(ltbmi)## model 1
					+ scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
					+ scale(lh_crp)  ##model 4
					+ strata(ID)
	,data)
	rst = rbind(rst, summary(model)$coef[1,])
}
rownames(rst) = S4_valid_measures

test = sample(which(data$inz_mi==1), 20)
train = data[-test,]
tr = train$inz_mi
x = train[,c("ltalteru", "lcsex", "my.diab","my.cigreg","my.alkkon")]
rst.match = Match(Tr = tr, X =x, exact = T, M = 4, ties = F)
train = train[c(rst.match$index.control, unique(rst.match$index.treated)), ]
#data$ID = rep(1:table(data$inz_mi)[1], 2)
train$ID = c(rep(1:table(train$inz_mi)[2], each = 4),1:table(train$inz_mi)[2])
require(survival)
rst = NULL
for(m in S4_valid_measures){
	train$m = log(train$Arg)-log(train$Trp)
	#train$m = log(train[,m])
	model = clogit(inz_mi ~ m 
					+ scale(ltbmi)## model 1
					+ scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
					+ scale(lh_crp)  ##model 4
					+ strata(ID)
			,train)
	rst = rbind(rst, summary(model)$coef[1,])
}
rownames(rst) = S4_valid_measures

candidate = rownames(rst)[which(rst[,5]<0.05)]
train$Arg_Trp = train$Arg/train$Trp
model = clogit(inz_mi ~ Ala + Arg_Trp +Gly+Met+Phe+Ser+Trp+lysoPC_a_C16_0+lysoPC_a_C17_0+lysoPC_a_C18_0+lysoPC_a_C18_2+PC_aa_C28_1+PC_ae_C36_4+PC_ae_C36_5+PC_ae_C38_5
				+ scale(ltbmi)## model 1
				+ scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
				+ scale(lh_crp)  ##model 4
				+ strata(ID)
		,train)

stepAIC(model)

model = clogit(inz_mi ~ Ala + Arg_Trp +Gly+lysoPC_a_C17_0+PC_ae_C36_5
				+ scale(ltbmi)## model 1
				+ scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
				+ scale(lh_crp)  ##model 4
				+ strata(ID)
		,train)

testset = data[c(test,which(data$inz_mi==0)),]
tr = testset$inz_mi
x = testset[,c("ltalteru", "lcsex", "my.diab","my.cigreg","my.alkkon")]
rst.match = Match(Tr = tr, X =x, exact = T, M = 1, ties = F)
testset = testset[c(rst.match$index.control, unique(rst.match$index.treated)), ]
testset$ID = c(rep(1:table(testset$inz_mi)[2], each = 1),1:table(testset$inz_mi)[2])

pred = predict(model, testset, type = "risk")
require(pROC)
plot(roc( testset$inz_mi, pred))