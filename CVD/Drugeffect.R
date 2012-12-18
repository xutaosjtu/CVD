# TODO: investigate the influence of statine on the candiate metabolites
# 
# Author: tao.xu
###############################################################################

data = data.frame(
		#time = S4$mi_time, event = S4$inz_mi,  # time and events
		ltmstati = as.factor(2-S4$ltmstati),
		scale(S4$ltalteru), scale(S4$ltbmi), as.factor(S4$lcsex), ##model 1
		S4$my.diab, ##model 2
		scale(S4$ltsysmm),  S4$my.cigreg, S4$my.alkkon, scale(S4$ll_hdla), scale(S4$ll_chola), ##model 3
		scale(log(S4$lh_crp)), ##model 4
		scale(log(S4[, metabo.asso])), scale(S4[,metabo.ratio.asso])
)
clinical = c("age", "ltbmi", "sex", "diabetes", "ltsysmm", "smoking", "alkkon", "ll_hdla", "ll_chola", "lh_crp")
colnames(data)[2:11] = clinical#, "total2HDL"
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

#############	effect of statine	################
rst=NULL
for(i in 1:length(metabo.ratio.asso)){
	model = glm(ltmstati ~. ,
	family = binomial(link = "logit")
	, data = data[, c(metabo.ratio.asso[i], clinical, "ltmstati")])
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rownames(rst) = metabo.ratio.asso


rst=NULL; candidates = c(metabo.asso, metabo.ratio.asso)
for(i in 1:length(candidates)){
	m = data[, candidates[i]]
	model = lm( 
			m ~ .,
			#family = binomial(link = "logit")
			subset = which(S4$inz_mi==1),
			data = data[, c("ltmstati", clinical[-10])]
	)
	
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rst = data.frame(rst, FDR = p.adjust(rst[,4], method = "BH"), bonferroni = p.adjust(rst[,4], method = "bonferroni"))
rownames(rst) = candidates
write.csv(rst, file = "Statine effect on metabolites and ratios.csv")

#############	association with c-reactive protein	#############
rst = NULL; candidates = c(metabo.asso, metabo.ratio.asso)
for(i in 1:length(candidates)){
	model = lm(lh_crp ~ ., data = data[, c(candidates[i], clinical)])
	rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = candidates


model = lm(lh_crp ~ ., data = data[, c(metabo.selected3, clinical)])
write.csv(summary(model)$coef[2:4,], file = "metabolite association with C reactive protein.csv")

















