require(metafor)

dir("meta analysis")

S2 = read.csv("meta analysis/S2_model1.csv",row.names=1)
S4 = read.csv("meta analysis/S4_model1.csv",row.names=1)

metabolites = intersect(rownames(S2), rownames(S4))
rst = NULL
for(i in metabolites){
  tmp = rbind(S2[i,], S4[i,1:5])
  tmp$ni = c(661,1342)
  rma.test = rma(yi = tmp$coef, sei= tmp$se.coef., measure="RR")
  summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
  rst = rbind(rst, summary)
}
rownames(rst) = metabolites
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst,"meta analysis/meta_fix_model1.csv")