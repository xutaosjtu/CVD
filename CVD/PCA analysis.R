# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


data = S4[which(S4$lthyact==2),c(valid_measures,S4.feature)]
adj = c("ltalteru", "lcsex", "ltbmi", "my.cigreg", "my.alkkon", "my.diab", "ll_chola", "ll_hdla")
residue(data, valid_measures, adj, control_group = which(data$lthya))

prcomp.hyper=prcomp( log(S4[which(S4$uthyact==1&S4$lthyact==2),c(valid_measures)]))

