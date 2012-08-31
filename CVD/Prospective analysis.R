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
for (m in S4_valid_measures){
	data = data.frame(surviv = S4[as.character(Cohort$zz_nr_s4), "CVD"], S4[as.character(Cohort$zz_nr_s4), m], S4[as.character(Cohort$zz_nr_s4),feature.cont], as.factor(S4[as.character(Cohort$zz_nr_s4),feature.disc]), t_time = apply(F4[as.character(Cohort$zz_nr_f4), c("utmialt","utschalt", "utalteru")], 1, min, na.rm = T)-S4[as.character(Cohort$zz_nr_s4),"ltalteru"])
	
	CAD.cox = coxph(Surv(t_time, surviv) ~  ., subset = index , data)
	
	
}
