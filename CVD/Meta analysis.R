require(metafor)

dir("meta analysis")

S2 = read.csv("meta analysis/S2_model4.csv",row.names=1)
S4 = read.csv("meta analysis/S4_model4.csv",row.names=1)

metabolites = intersect(rownames(S2), rownames(S4))
#pdf("meta analysis/meta_fixed_model4.pdf")
par(mfrow=c(3,3))
rst = NULL
for(i in metabolites){
  tmp = rbind( S4[i,1:5], S2[i,])
  tmp$ni = c(1342,661)
  rma.test = rma(yi = tmp$coef, sei= tmp$se.coef., 
                 method = "FE", 
                 measure="RR")
  #forest(rma.test, main = i)
  summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
  rst = rbind(rst, summary)
}
#dev.off()

rownames(rst) = metabolites
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
rst[,"fdr"] = p.adjust(rst[,"pvalue"], method = "BH")
write.csv(rst,"meta analysis/meta_fix_model4.csv")

## Meta-analysis of CRP effect size in S4 and S2
S2 = read.csv("estimates of confounders in S2_model 4.csv", row.names = 1)
S4 = read.csv("estimates of confounders in S4_model 4.csv", row.names = 1)

rst = NULL
for(i in 1:nrow(S2)){
  tmp = rbind(S2[i,], S4[i,])
  rma.test = rma(yi = tmp$coef, sei = tmp$se.coef., 
                 #method = "FE",
                 measure = "RR")
  summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
  rst = rbind(rst, summary)
}
rownames(rst) = rownames(S2)
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, "meta analysis_random_confounders model4.csv")  
