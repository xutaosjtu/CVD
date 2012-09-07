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


###	new developped CVDs
index = which(F4[as.character(Cohort$zz_nr_f4),"CVD"])
data.newCVD = data.frame(S4[as.character(Cohort$zz_nr_s4)[index], c("CVD", "ltalteru") ] , F4[as.character(Cohort$zz_nr_f4)[index],c("CVD","utmialt","utschalt","utalter")])
write.csv(data.newCVD, file = "newly developped CVDs.csv")

index2 = which(S4[as.character(Cohort$zz_nr_s4),"CVD"])
data.newCVD = data.frame(S4[as.character(Cohort$zz_nr_s4)[index2], c("CVD", "ltalteru") ] , F4[as.character(Cohort$zz_nr_f4)[index2],c("CVD","utmialt","utschalt","utalter")])
write.csv(data.newCVD, file = "CVDs in S4.csv")

table(F4[as.character(Cohort$zz_nr_f4)[-c(index, index2)],"CVD"])
table(S4[as.character(Cohort$zz_nr_s4)[-c(index, index2)],"CVD"])

# assessing the biomarkers for prediction of CVD (on 26 cases!!!) 
index = which(!(S4[as.character(Cohort$zz_nr_s4),"CVD"]))
substr(feature.cont, 1, 2) <- "l"
substr(feature.disc, 1, 2) <- "l"
rst = NULL;
for (m in S4_valid_measures){
	data = data.frame(
			CVD = F4[as.character(Cohort$zz_nr_f4), "CVD"],
			
			CVDS4 = S4[as.character(Cohort$zz_nr_s4), "CVD"], 
			
			S4[as.character(Cohort$zz_nr_s4), m],
			
			S4[as.character(Cohort$zz_nr_s4),feature.cont],
			
			S4[as.character(Cohort$zz_nr_s4),feature.disc], 
			
			t_time = apply(F4[as.character(Cohort$zz_nr_f4), c("utmialt","utschalt", "utalteru")], 1, min, na.rm = T)-S4[as.character(Cohort$zz_nr_s4),"ltalteru"]
	)
	
	data = na.omit(data)
	
	index = which(!data$CVDS4 & data$t_time>0)
	
	colnames(data) = c("CVD" , "CVDS4" , "metabolite", "ltalteru" , "ltsysmm" , "ll_hdln" , "ll_choln" , "lp_diab_who06" , "lcsex" , "t_time")
	
	CVD.cox = coxph(Surv(t_time, CVD) ~  log(metabolite) + ltalteru + log(ltsysmm) + log(ll_hdln) + log(ll_choln) + as.factor(lp_diab_who06) + as.factor(lcsex) , subset = index , data)
	
	rst = rbind(rst , summary(CVD.cox)$coefficients[1,])
}
rownames(rst) = S4_valid_measures


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

Comparison.prosective<- function(baseline, feature, metabo, adj, subset){
#	baseline	---		data frame of variables at baseline 
#	feature		---		list or matrix indicating the disease state in at different time points
#	metabo		---		names of metabolites
#
	for(i in 1:length(metabo)){
		
		data = data.frame(log(baseline[, metabo[i]]),
				basline[,adj], feature)
		
		model = glm(interaction(feature[,1], feature[,2]) ~ . , subset, data, family = binomial(link = "logit"))
		
		rst = rbind(rst, summary(model)$coefficients[2, ])
	
	}
	
	rownames(rst) = metabo
	
	return (rst)

}



for(i in 1:4){
	
	tmps4 = S4[Cohort$zz_nr_s4, diseases.s4[i]]
	tmpf4 = F4[Cohort$zz_nr_f4, diseases.f4[i]]
	feature = data.frame(tmps4, tmpf4)
	subset = intersect (which(feature[,1] == 2), which(feature[,2] !=3 ))
	rst = Comparison.prospective(S4[], feature, metabo = S4_valid_measures, subset = subset, adj = model_1)
	
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
}












