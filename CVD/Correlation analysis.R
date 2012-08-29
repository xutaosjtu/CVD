# Correlation analysis of metabolites and other clinical parameters
# 
# Author: tao.xu
###############################################################################

# pearson correlation
clinical = scan(what = character())
ll_choln
ll_hdln
ll_ldln
ll_hgb
ll_hbav
ll_trin
lh_crp

substr(clinical, 1, 2) <- "l"
data = data.frame(S4[,clinical], S4[,S4_valid_measures])
data <- na.omit(data)
#data = log(data)
cor.pearson = cor(data, use = "pairwise.complete.obs") 

# partial correlation
require(corpcor)
cor.par = cor.partial (data, colnames(data))
diag(cor.par) = 0
# partial correlation with shrinkage
cor.shrink = pcor.shrink(data)

require(gplots)
pdf("Correlation of metabolites.pdf", width = 20, height = 20)
heatmap.2(cor.pearson, col = greenred(200)[60:200], trace = "none", keysize  = 1, main = "Pearson correlation")
heatmap.2(cor.par, col = greenred(200)[32:193], trace = "none", keysize  = 1, main = "Partial correlation")
heatmap.2(cor.shrink, col = greenred(200)[40:200], trace = "none", keysize  = 1, main = "Partial correlation with shrinkage")
dev.off()

net = cor2link(Zeta = cor.shrink, threshold = 0.12)
write.csv(net, "association of metabolites with clinical measurements_0.12_hba1c2.csv", row.names = F)

######################	F4	######################
substr(clinical, 1, 2) <- "u"
F4_valid_measures = gsub(".", replacement= "_", F4_valid_measures, fixed = T)

data = data.frame(F4[,clinical], F4[,F4_valid_measures])
data <- na.omit(data)
cor.pearson = cor(data, use = "pairwise.complete.obs") 

# partial correlation
require(corpcor)
cor.par = cor.partial (data, colnames(data))

# partial correlation with shrinkage
cor.shrink = pcor.shrink(data)

require(gplots)
pdf("Correlation of metabolites F4.pdf", width = 20, height = 20)
heatmap.2(cor.pearson, col = greenred(200)[60:200], trace = "none", keysize  = 1, main = "Pearson correlation")
heatmap.2(cor.par, col = greenred(200)[32:193], trace = "none", keysize  = 1, main = "Partial correlation")
heatmap.2(cor.shrink, col = greenred(200)[40:200], trace = "none", keysize  = 1, main = "Partial correlation with shrinkage")
dev.off()

net = cor2link(Zeta = cor.par, threshold = 0.13)
write.csv(net, "F4 association of metabolites with clinical measurements_0.13_hba1c.csv", row.names = F)
