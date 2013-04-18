# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


data = S4[which(S4$lthyact==2),c(valid_measures,S4.feature,"lthyact","uthyact")]
adj = c("ltalteru", "lcsex", "ltbmi", "my.cigreg", "my.alkkon", "my.diab", "ll_chola", "ll_hdla")
data.resid=residue(data, valid_measures, adj, control_group = which(data$lthyact==2&data$uthyact==2))

data.resid = preprocess(data.resid, Metabolites = colnames(data.resid))
prcomp.hyper=prcomp(data.resid[which(data$lthyact==2&data$uthyact==1),])

tmp = as.matrix(data.resid[which(data$lthyact==2&data$uthyact==1),])
heatmap(tmp)
dim(tmp)
tmp = apply(tmp,2,scale)