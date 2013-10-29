##########################
##  expression analysis ##
##########################
## association between lpcs and LPPLA2 gene
subset = which(!is.na(F4$zz_nr_f4_bio) & F4$expr_in_F4!='')
subset.exp = F4$zz_nr_s4f4_genexp[subset]
subset.metab  = F4$zz_nr_f4_bio[subset]
F4.sub = F4[subset, ]
F4.sub$expression = F4.expression["ILMN_1701195", as.character(subset.exp)]
#metab = F4.metab[as.character(subset.metab), "lysoPC.a.C17.0"]
rst = NULL
for(i in c(125:137)){
  F4.sub$metab = F4.metab[as.character(subset.metab), i]
  model = lm(expression ~ metab + utalteru + as.factor(ucsex) + as.factor(my.diab), data = F4.sub)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = colnames(F4.metab)[c(125:137)]

#S4
subset = which(!is.na(S4$zz_nr_s4_bio) & S4$expr_in_S4!="")
subset.exp = S4$zz_nr_s4f4_genexp[subset]
subset.metab  = S4$zz_nr_s4_bio[subset]
S4.sub = S4[subset,]
S4.sub$expression = S4.expression["ILMN_1701195", as.character(subset.exp)]
rst = NULL
for(i in c(105:118)){
  S4.sub$metab = S4.metab[as.character(subset.metab), i]
  model = lm(expression ~ scale(metab) + ltalteru + as.factor(lcsex), data = S4.sub)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = colnames(S4.metab)[c(105:118)]

# F3
subset = which(!is.na(F3$zz_nr_f3_genexp))
F3.sub = F3[subset,]
F3.sub$expression = t(F3.expression["ILMN_1701195",as.character(F3.sub$zz_nr_f3_genexp)])
rst = NULL
for(i in c(175:189)){
  F3.sub$metab = F3.sub[, i]
  model = lm(expression ~ scale(metab) + rtalteru + as.factor(rcsex), data = F3.sub)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = colnames(F3)[c(175:189)]


S4.F4 = read.csv("metabolites/K2912_Xu_S4F4_trans120213.csv", sep=";")
F4 = S4.F4[, c(34:70)]
F4.metab = read.csv(file="metabolites/KoraF4-metabolomics-quality.controlled-mice.imputed-20100107.csv",sep=";", row.names=1)
rownames(F4)=F4$sample.id

load("expression/Expression_f4_adjusted_technical_variables.Rdata")
F4.expression = adj.expr.f4
rm(adj.expr.f4)

subset = which(!is.na(F4$zz_nr_f4_bio) & F4$expr_in_F4!='')
subset.exp = F4$zz_nr_s4f4_genexp[subset]
subset.metab  = F4$zz_nr_f4_bio[subset]
F4.sub = F4[subset, ]

require(doMC)
registerDoMC(cores = 8)

association = foreach(i = 1:nrow(F4.expression), .combine = rbind ) %dopar% {
  F4.sub$expression = F4.expression[i, as.character(subset.exp)]
  rst.tmp = NULL
  for(j in 124:136){
    F4.sub$metab = F4.metab[as.character(subset.metab), j]
    model = lm(expression ~ metab + utalteru + as.factor(ucsex), data = F4.sub)
    rst.tmp = c(rst.tmp, summary(model)$coef[2,]) 
  }
  return(rst.tmp)
}

load("LPCs metabolite associated genes.RData")
ls()
names(metab_associated_gene) = colnames(F4.metab)[124:136]

lapply(metab_associated_gene, 
       function(x){
         return(annotation[which(annotation[,1] %in% names(x)),2])
       }
)
