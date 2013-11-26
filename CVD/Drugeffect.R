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
		scale(log(S4[, metabo.asso])), scale(S4[,metabo.ratio.asso]),
		prev_mi=S4$prev_mi
)
clinical = c("ltalteru", "ltbmi", "lcsex", "my.diab", "ltsysmm", "my.cigreg", "my.alkkon", "ll_hdla", "ll_chola", "lh_crp")
colnames(data)[2:11] = clinical#, "total2HDL"
na.index = unique(unlist(apply(data, 2, function(x) which(is.na(x)))))

#############	effect of statine	################
rst=NULL; candidates = c(metabo.asso, metabo.ratio.asso)
for(i in 1:length(candidates)){
	model = glm(ltmstati ~. ,
	family = binomial(link = "logit")
	, data = data[, c(candidates[i], clinical, "ltmstati")])
	rst = rbind(rst, summary(model)$coefficients[2,])
}
rownames(rst) = candidates


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
rst = NULL; 
for(i in metabo.asso){
  S4$m = scale((S4[, i]))
	model = lm(log(lh_crp) ~ m + 
               scale(ltalteru) + as.factor(lcsex)
	           + scale(ltbmi) + as.factor(my.diab)  ##model 2
	           + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3
	           , subset = which(S4$prev_mi==0) 
             , data = S4)
	rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = metabo.asso
rst = data.frame(rst, lower = rst[,1]-1.96*rst[,2], upper = rst[,1]+1.96*rst[,2])
write.csv(rst, file = "Associations of candidates with CRP_S4.csv")


rst = NULL; 
for(i in metabo.asso){
  S2$m = scale(log(S2[, i]))
  model = lm(log(cl_crp) ~ m + 
               scale(ctalteru) + as.factor(ccsex)
             + scale(ctbmi)## model 1
             + as.factor(my.diab)  ##model 2
             + scale(ctsysmm) + as.factor(my.cigreg) + my.alkkon  + scale(cl_chola) + scale(cl_hdla)
             , subset = prev_mi02==0&subcoho==1| inz_mi02==1
             , data = S2)
  rst = rbind(rst, summary(model)$coef[2,])
}
rownames(rst) = metabo.asso
rst = data.frame(rst, lower = rst[,1]-1.96*rst[,2], upper = rst[,1]+1.96*rst[,2])
write.csv(rst, file = "Associations of candidates with CRP_S2.csv")


model = lm(scale(log(lh_crp)) ~ #scale(log(Arg)) + scale(log(Trp)) + scale(log(lysoPC_a_C17_0)) + scale(log(PC_aa_C32_2)) + scale(log(lysoPC_a_C18_2)) + 
            scale(ltalteru) + as.factor(lcsex)
           + scale(ltbmi) + as.factor(my.diab)  ##model 2
           + scale(ltsysmm) + as.factor(my.cigreg) + as.factor(my.alkkon)  + scale(ll_chola) + scale(ll_hdla) ##model 3
           ,subset = which(S4$prev_mi==0),
           data = S4)


plot(lysoPC_a_C17_0 ~ log(lh_crp), data = S4, subset = which(S4$prev_mi==0))

require(grid)

model = lm(lh_crp ~ ., subset = which(S4$prev_mi==0), data = S4[, c(clinical)])
count=0;rcount=0
for(i in metabo.selected3){
	#data = data.frame(marker = S4[,i],crp = log(S4$lh_crp))
	indx = as.vector(model$na.action)
	#tmp =data.frame(marker = S4[which(S4$prev_mi==0)[-indx],i], crp = model$residuals)
	tmp = data.frame(marker = S4[which(S4$prev_mi==0),i], crp = S4[which(S4$prev_mi==0),"lh_crp"])
	x_min = quantile(tmp$marker,0.01)
	x_max = quantile(tmp$marker,0.99)
	
	tmp = tmp[which(tmp$marker<x_max&tmp$marker>x_min),]
	
	c <- ggplot(tmp ,aes(marker,crp))+theme_bw()
	c <- c  + geom_point(alpha = 0.5) + xlab(i) +ylab("log(CRP)")+ stat_smooth(size =1,color="red",method='lm')
	
	###make different panels for plots+ coord_cartesian(xlim = c(x_min,x_max))
	pushViewport(
			viewport(x=0.25+count*5/10, y=0.75-(rcount*5)/10,width = 0.5, height = 0.5, angle = 0)
	)
	print(c,newpage=F)
	upViewport() 
	
	####format control, print 3 figures for each row
	count=count+1;
	if(count%%2==0){rcount=rcount+1;count=0}
}












