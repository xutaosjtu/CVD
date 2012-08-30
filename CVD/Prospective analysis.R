# TODO: Add comment
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
