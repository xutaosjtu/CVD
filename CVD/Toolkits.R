# Functions could be used in the analysis
# 
# Author: tao.xu
###############################################################################

#function(data, subset = NULL, feature, disease, files){
#	if(is.null(subset)){
#		subset = dim(S4)[1]
#	}
#	for(i in 1:length(diseases)){
#		chars = apply(S4[subset,feature], 2, function(x) tapply(x, INDEX = as.factor(S4[subset, diseases[i]]), mean, na.rm= T))
#		write.csv(t(chars), file = paste(names(diseases)[i],".csv", sep = ""))
#	}
#}
#testfun = function(x, index, FUN){
# 	tapply(x, INDEX = index, FUN)
#}


characteristics = function(data , factor, d, ...)
{
	chars = apply(
			data, 2, 
			function(x) {
				m = tapply(x, INDEX = factor, mean, ...)
				s = tapply(x, INDEX = factor, sd, ...)
				m = round(m , digits = d)
				s = round(s, digits = d)
				return(paste(m, " (", s, ")", sep = ""))
			}
	)
	return(chars)
}


# partial correlation calculation
cor.partial <- function(data, variables, ...)
{
	P = cor( data[,variables], ... )
	P.invers = solve( P )
	P.diag = diag( P.invers )
	Z = - P.invers / sqrt(P.diag %*% t(P.diag))
	return(Z)
}

#correlation matrix to linkage pairs
cor2link = function(Zeta, threshold)
{
	links = NULL;
	for(i in 1:(dim(Zeta)[1]-1)){
		for(j in (i+1):dim(Zeta)[1]){
			if(abs( Zeta[i,j] ) >= threshold){
				tmp = c( colnames( Zeta )[i], colnames( Zeta )[j], 
						Zeta[i,j])
				links = rbind(links, tmp)
			}
		}
	}
	return(links)
}

#	Differential correxpression
#	The method was described in Cho, Sung, Jihun Kim, and Ju Kim. “Identifying Set-wise Differential Co-expression in Gene Expression Microarray Data.” BMC Bioinformatics 10, no. 1 (2009): 109.
#	Description of parameters:
#	cor1, cor2: two pearson correlation values
#	N1, N2: the number of samples used for the calculation of cor1 and cor2
diffcorr <- function(cor1, cor2, N1, N2)
{
	Z1 = 0.5*log((1+cor1)/(1-cor1))
	Z2 = 0.5*log((1+cor2)/(1-cor2))
	p = 1 - pnorm(abs((Z1-Z2) / sqrt(1/(N1-3) + 1/(N2-3))))
	return(p)
}

## coversion of parameters between S4 and F4
#function(parameters , optional = "l"){
#	substr(parameters, 1, 2) <- optional
#}

# regression analysis
logisticRegression = function(meta , disease, valid_measures , feature.cont, feature.disc, metalog = TRUE, ...)
{
	rst = NULL 
	fdisc = NULL; fcont = NULL
	for (f in feature.disc){
		fdisc = cbind(fdisc, as.factor(meta[,f]))
	}
	for(m in valid_measures){
		if(metalog) { 
			metabolite = log(meta[,m])
			fcont = log(meta[ , feature.cont])
		}
		else { 
			metabolite = meta[,m]
			fcont = meta[ , feature.cont]
		}
		data = data.frame (disease, metabolite , fcont, fdisc)
		data = na.omit(data)
		model = glm(as.factor(disease) ~ ., data = data, family = binomial(link = "logit"), na.action = na.omit , ...)
		rst = rbind(rst , summary(model)$coefficients[2,])
	}
	rst = data.frame(rst, p.adjust(rst[,4], method = "BH"))
	rownames(rst) = valid_measures
	return(rst)
}

#########	longitudinal analysis using logistic regression	################
Comparison.prospective<- function(baseline, feature, metabo, adj, subset){
#	baseline	---		data frame of variables at baseline 
#	feature		---		list or matrix indicating the disease state in at different time points
#	metabo		---		names of metabolites
#
	rst = NULL
	
	for(i in 1:length(metabo)){
		
		data = data.frame(log(baseline[, metabo[i]]), baseline[,adj])
		
		model = glm(interaction(feature[,1], feature[,2]) ~ . , data = data, subset = subset,  family = binomial(link = "logit"))
		
		rst = rbind(rst, summary(model)$coefficients[2, ])
		
	}
	print(dim(rst)); print(i)
	rownames(rst) = metabo
	
	return (rst)
	
}

