# regression analysis
#	model 1: basic model
#	feature.cont = c("ltalteru", "ltalkkon" , "ltbmi")
# 	feature.disc = c("lcsex")
#	model 2: full model
#	feature.cont = 
#	feature.disc = c("lcsex")
#	model 3: Estimated clinical model
#	Hypertension:
#	feature.cont = c("ltalteru","ll_hgb", "ll_trin", "ltbmi" , "lcsex" , "lh_crp", "ll_hdln") 
#	feature.disc = c("lcsex")
#	Stroke: 
#	feature.cont = c("ltalteru" , "lh_crp" , "ll_choln" )
#	feature.disc = c("lcsex")
#	Myocardial infarction: 
#	feature.cont = c("ll_choln" , "ll_ldln" , "lthumf" )
#	feature.disc = c("lcsex")
#	Heart failure
#	feature.cont = c("ltalteru"	, "ll_choln" , "ll_hdln" , "ll_trin" , "lthumf" , "lttumf" , "ltalkkon" , "ll_hgb")
#	feature.disc = c("lcsex")
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
	data = data.frame (disease = as.factor(S4[,diseases[1]]), metabolite = log(S4[,m]), S4[,feature.cont], as.factor(S4[,feature.disc]))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2:4,4] )
}
rst = data.frame(rst, dis1_FDR = p.adjust(rst[,1], method = "BH"), dis2_FDR = p.adjust(rst[,2], method = "BH"), dis3_FDR = p.adjust(rst[,3], method = "BH"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = paste( names(diseases)[1], "categorical model_3 S4.csv", sep = " "))

rst = NULL
for (m in S4_valid_measures){
	data = data.frame (disease = S4[,diseases[1]], metabolite = log(S4[,m]), S4[,feature.cont], as.factor(S4[,feature.disc]))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2,] )
}
rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
rownames(rst) = S4_valid_measures
write.csv(rst, file = paste( names(diseases)[1], "model_3 S4.csv", sep = " "))


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
	data = data.frame (disease = as.factor(F4[,diseases[1]]), metabolite = log(F4[,m]), F4[,feature.cont], as.factor(F4[,feature.disc]))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2:4,4] )
}
rst = data.frame(rst, dis1_FDR = p.adjust(rst[,1], method = "BH"), dis2_FDR = p.adjust(rst[,2], method = "BH"), dis3_FDR = p.adjust(rst[,3], method = "BH"))
rownames(rst) = F4_valid_measures
write.csv(rst, file = paste( names(diseases)[1], "categorical model_3 F4.csv", sep = " "))

rst = NULL
for (m in F4_valid_measures){
	data = data.frame (disease = F4[,diseases[1]], metabolite = log(F4[,m]), F4[,feature.cont], as.factor(F4[,feature.disc]))
	model = lm( metabolite ~ ., data = data)
	rst = rbind ( rst , summary(model)$coefficients[2,] )
}
rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
rownames(rst) = F4_valid_measures
write.csv(rst, file = paste( names(diseases)[1], "model_3 F4.csv", sep = " "))

# stroke, myocardio infarction, heart failure
for(i in 4:4){
	rst = NULL
	for(m in F4_valid_measures){
		data = data.frame (disease = F4[,diseases[i]], metabolite = log(F4[,m]), F4[,feature.cont], as.factor(F4[,feature.disc]))
		
		model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), subset = which(data$disease!=3))
		rst = rbind ( rst , summary(model)$coefficients[2,] )
	}
	rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
	rownames(rst) = F4_valid_measures
	write.csv(rst, file = paste( names(diseases)[i], "model_3 F4.csv", sep = " " ))
}	


#################
tmp = (S4[,diseases[2:4]])
D = apply(tmp, 1, function(x) sum(x==1))
table(D)
S4$CVD =  apply(tmp, 1, function(x) sum(x==1))!=0

tmp = (F4[,diseases[2:4]])
table(apply(tmp, 1, function(x) sum(x==1)))


