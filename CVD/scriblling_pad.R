# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


time <-read.csv("F:/Cardiovascular disease/data/sas data/pv95_12_add.csv")


setwd("H:/Met_ED_Val/")
load("data.RData")
require(boot)

var_f = function(x, inds){sd(x[inds], na.rm = T)/mean(x[inds], na.rm = T)}

indx = which(!duplicated(Met_ED_Val[,c("bmi_dat","bmi")]))
indx1 = which(Met_ED_Val$X!=2&Met_ED_Val$X!=3)
indx2 = which(Met_ED_Val$X==2|Met_ED_Val$X==3)

rst.edval = apply(Met_ED_Val[indx1,c(179,16:178)], 2, function(x) boot(x, var_f, R=100))
rst.edval2 = apply(Met_ED_Val[indx2,c(179,16:178)], 2, function(x) boot(x, var_f, R=100))

#indx = which(!duplicated(Met_Val[,c("bmi_dat","bmi")]))
#rst.val = apply(Met_Val[indx,15:177], 2, function(x) boot(x, var_f, R=100))

pdf("D:/CVs in different metabolites.pdf", width = 9, height = 10)
par(mfrow = c(3,5))
for( i in names(rst.edval)){
	boxplot(data.frame('fast' = rst.edval[[i]]$t, 'non-fast' = rst.edval2[[i]]$t), main = i, ylab = "Coefficient of variation")
}
dev.off()

### individual differences
mean.nfast = Met_ED_Val[which(Met_ED_Val$X==1)[1:46],c(16:178)] + Met_ED_Val[which(Met_ED_Val$X==4),c(16:178)] + Met_ED_Val[which(Met_ED_Val$X==5),c(16:178)]
mean.nfast = rbind(mean.nfast, Met_ED_Val[which(Met_ED_Val$X==1)[47:50]*3,c(16:178)])
mean.nfast = mean.nfast/3

mean.fast = Met_ED_Val[which(Met_ED_Val$X==2),c(16:178)] + Met_ED_Val[which(Met_ED_Val$X==3),c(16:178)]
mean.fast = mean.fast/2

rst = mean.nfast-mean.fast
pdf("D:/distribution of concentration differences in fasting and nonfasting samples.pdf", height = 12.5, width = 10)
par(mfrow = c(5,4))
for(i in colnames(rst)){
	hist(rst[,i], main = i, xlab = "non-fasting concentration- fasting concentration")
	abline(v = 0, col = "red", lwd = 2)
}
dev.off()

pdf("concentration of fasting against non-fasting.pdf", height = 18, width = 14)
par(mfrow = c(5,4))
for(i in colnames(rst)){
	pcc = cor(mean.nfast[,i],mean.fast[,i], use= "p", method ="s")
	if(i %in% F4_valid_measures){flag = ""}
	else {flag = "X"}
	plot(mean.nfast[,i], mean.fast[,i], main = paste(flag, i, "| spearman's rho=",round(pcc, 3), collapse=" "), ylab = "non-fasting", xlab="fasting")
	abline(lm(nfast~fast, data = data.frame('nfast' = mean.fast[,i], 'fast' = mean.nfast[,i])), col = "red")
	abline(0,1,lty = 2)
	legend("bottomright", c("if equal","linear regression"),  lty=c(2,1), col = c("black","red"))
	#lines(lowess(mean.nfast[,i]~ mean.fast[,i]))
}
dev.off()

pcc=NULL
for(i in names(rst)){
	pcc = c(pcc, cor(mean.nfast[,i],mean.fast[,i], use= "p", method ="s"))
}
