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


######## correlation between candidates with CRp
metabo.asso = scan(what = character())
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
PC_ae_C36_1
PC_ae_C36_2
PC_ae_C38_0
PC_ae_C38_2
PC_ae_C40_1
SM_C24_1

##GGM 
require(corpcor)
require(corrplot)
pdf("correlation with conventional biomarkers_S4.pdf")
tmp = S4[, c(metabo.asso, "lh_crp", "ltbmi",   "ltsysmm","ll_ldla", "ll_hdla")]#"ll_hbav",
#pdf("correlation with conventional biomarkers_S2.pdf")
#tmp = data.S2[, c(metabo.asso, "cl_crp", "ctbmi",  "ctsysmm","cl_ldla", "cl_hdla")]
tmp = log(tmp)
tmp = cor(tmp, use = "p" )
#tmp = cor2pcor(tmp)
rownames(tmp)=c(metabo.asso, "lh_crp", "ltbmi", "ltsysmm","ll_ldla", "ll_hdla")
#rownames(tmp) = c(metabo.asso,"cl_crp", "ctbmi",  "ctsysmm","cl_ldla", "cl_hdla")
colnames(tmp)=rownames(tmp)
diag(tmp) = 0
col3 <- colorRampPalette(c("blue", "white", "red"))
corrplot(tmp, col = col3(200), tl.col="black",
          type = "lower", tl.pos="l")
dev.off()
##partial correaltion adjusting for age and sex
# tmp = S4[, c(metabo.asso, "lh_crp", "ltbmi",  "ll_hbav", "ltsysmm","ll_ldla", "ll_hdla", "ltalteru","lcsex")]
# tmp = na.omit(tmp)
# tmp.metab = log(tmp[, c(metabo.asso)])
# tmp.biom = log(tmp[, c("lh_crp", "ltbmi",  "ll_hbav", "ltsysmm","ll_ldla", "ll_hdla")])
# tmp.adj = tmp[, c("ltalteru","lcsex")]
# tmp.adj$lcsex = 2-tmp.adj$lcsex

tmp = data.S2[, c(metabo.asso, "cl_crp", "ctbmi", "ctsysmm","cl_ldla", "cl_hdla", "ctalteru","ccsex")]
tmp = na.omit(tmp)
tmp.metab = log(tmp[, c(metabo.asso)])
tmp.biom = log(tmp[, c("cl_crp", "ctbmi", "ctsysmm","cl_ldla", "cl_hdla")])
tmp.adj = tmp[, c("ctalteru","ccsex")]
tmp.adj$ccsex = 2-tmp.adj$ccsex
tmp = sapply(tmp.biom, 
       function(x){
         resid.x = lm(x~., data = data.frame(x,tmp.adj))$residual
         rst = sapply(tmp.metab, 
                      function(y){
                        resid.y=lm(y ~., data = data.frame(y,tmp.adj))$residual
                        pcc =cor(resid.x,resid.y, use="pair")
                        #pcc = cor2pcor(pcc)
                        return(pcc)
                      }
                      )
         #rownames(rst) = c("cor", "upper", "lower", "p-value")
         return(rst)
       }
)
col3 <- colorRampPalette(c("blue", "white", "red"))
corrplot(tmp, col = col3(200), tl.col = "black")

require(deal)
tmp = S4[, c(metabo.selected, "lh_crp", "inz_mi")]
tmp$inz_mi=as.factor(tmp$inz_mi)
dim(tmp)
tmp = na.omit(tmp)
tmp.nw = network(tmp)
tmp.j = jointprior(tmp.nw, 12)

tmp.learn = learn(tmp.nw, tmp, tmp.j)
tmp.nw.learn = getnetwork(tmp.learn)

drawnetwork(tmp.nw, tmp, tmp.j)$nw
tmp.h<-heuristic(tmp.nw.learn, tmp, tmp.j,
                 trace = FALSE)

require(bnlearn)
tmp = S4[, c(metabo.selected, "lh_crp", "inz_mi")]
tmp = na.omit(tmp)
res = gs(tmp)
plot(res)
res2 = iamb(tmp)
plot(res2)

res = gs(tmp, debug=TRUE)
sink(file = "Casual inference/Bayesian network/learning.txt")