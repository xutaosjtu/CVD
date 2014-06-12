require(metafor)

dir("meta analysis")

S2 = read.csv("meta analysis/S2_model1.csv",row.names=1)
S4 = read.csv("meta analysis/S4_model1.csv",row.names=1)

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

refine = read.csv("meta analysis/REFINE_model1.csv", row.names = 1)
refine$se.coef. = (log(refine$upper)-log(refine$OR))/1.96
refine$coef = log(refine$OR)

require(metafor)
## meta analysis of KORA S4 and refine
rst = NULL
for(m in candidates){
  tmp = rbind(S4[m,c("coef", "se.coef.")], refine[m,c("coef", "se.coef.")])
  tmp = as.data.frame(tmp)
  rma.test = rma(yi = tmp$coef, sei = tmp$se.coef., n1i = 1342, n2i = 108*3,
                 #method = "FE",
                 measure = "RR")
  rst = rbind(rst , c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp))
}
rownames(rst) = candidates
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, file = "meta-analysis_model 3_refine_S4.csv")

## meta analysis of KORA S2 and refine
rst = NULL
for(m in candidates){
  tmp = rbind(S2[m,c("coef", "se.coef.")], refine[m,c("coef", "se.coef.")])
  tmp = as.data.frame(tmp)
  rma.test = rma(yi = tmp$coef, sei = tmp$se.coef., n1i = 1342, n2i = 108*3,
                 #method = "FE",
                 measure = "RR")
  rst = rbind(rst , c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp))
}
rownames(rst) = candidates
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, file = "meta-analysis_model 3_refine_S2.csv")

## meta analysis of KORA S4, S2 and refine
rst = NULL
for(m in candidates){
  tmp = rbind(S4[m,c("coef", "se.coef.")], S2[m,c("coef", "se.coef.")],refine[m,c("coef", "se.coef.")])
  tmp = as.data.frame(tmp)
  rma.test = rma(yi = tmp$coef, sei = tmp$se.coef., n1i = 1342, n2i = 108*3,
                 #method = "FE",
                 measure = "RR")
  pdf(paste(m,"_model 1_meta.pdf"))
  forest(rma.test, main = m)
  dev.off()
  rst = rbind(rst , c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp))
}
rownames(rst) = candidates
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, file = "meta-analysis_model 3_refine_S4_S2.csv")



## Meta-analysis of CRP effect size in S4 and S2
S2 = read.csv("CRP association analysis/Change in CRP effect size/estimates of confounders plus 5 metabolites in S2_model 4.csv", row.names = 1)
S4 = read.csv("CRP association analysis/Change in CRP effect size/estimates of confounders plus 5 metabolites in S4_model 4.csv", row.names = 1)
refine = read.csv("CRP association analysis/Change in CRP effect size/estimates of confounders plus 5 metabolites in REFINE_model 4_cl.csv", row.names = 1)
#S2.estimate = read.csv("Change in CRP effect/estimates of confounders plus original four metabolites in S2_model 4.csv", row.names = 1)
#S4.estimate = read.csv("Change in CRP effect/estimates of confounders plus original four metabolites in S4_model 4.csv", row.names = 1)

refine$se.coef. = (log(refine$Upper)-log(refine$HR))/1.96
refine$coef = log(refine$HR)

rst = NULL
for(i in 1:nrow(S2)){
  tmp = rbind(S4[i,c("coef", "se.coef.")]
              ,S2[i,c("coef", "se.coef.")]
              ,refine[i,c("coef", "se.coef.")]
              )
  tmp = as.data.frame(tmp)
  rma.test = rma(yi = tmp$coef, sei = tmp$se.coef., 
                 #method = "FE",
                 measure = "RR")
  summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
  rst = rbind(rst, summary)
}
rownames(rst) = rownames(S2)
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, "meta analysis_random_S4 S2 refine_metabo+model4_cl.csv")  

rma.baseline = rma.test 
tmp = rbind(S4.estimate [nrow(S4.estimate ),], S2.estimate [nrow(S2.estimate ),1:5])
rma.plusmetab = rma(yi = tmp$coef, sei = tmp$se.coef., 
                    #method = "FE",
                    measure = "RR")


tmp = argtrp[,1:5]
tmp$ni = c(1342,661)
rma.test = rma(yi = tmp$coef, sei= tmp$se.coef., 
               method = "FE", 
               measure="RR")
pdf("Arg Trp ratio model 3.pdf")
forest(rma.test)
dev.off()



## linear associations between the metabolites with CRP
S2 = read.csv("CRP association analysis/Associations of candidates with CRP_S2.csv", row.names = 1)
S4 = read.csv("CRP association analysis/Associations of candidates with CRP_S4.csv", row.names = 1)
refine = read.csv("CRP association analysis/Associations of candidates with CRP_Refine.csv", row.names = 1)

rst = NULL
for(i in 1:nrow(S2)){
  tmp = rbind(S4[i,c("Estimate", "Std..Error")]
              ,S2[i,c("Estimate", "Std..Error")]
              ,refine[i,c("Estimate", "Std..Error")]
  )
  tmp = as.data.frame(tmp)
  colnames(tmp) = c("yi", "sei")
  rma.test = rma(yi= tmp$yi, sei=tmp$sei, method = "FE")
  summary = c(rma.test$b, rma.test$se, rma.test$ci.lb, rma.test$ci.ub,rma.test$pval, rma.test$QE, rma.test$QEp)
  rst = rbind(rst, summary)
}
rownames(rst) = rownames(S2)
colnames(rst) = c("estimates", "se", "ci.lb","ci.ub", "pvalue","heterogenity", "Hetero Pvalue")
write.csv(rst, "CRP association analysis/Associaitons of candidates with CRP_meta analysis_random_S4 S2 refine.csv")  
