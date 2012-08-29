# Functions could be used in the analysis
# 
# Author: tao.xu
###############################################################################

function(data, subset = NULL, feature, disease, files){
	if(is.null(subset)){
		subset = dim(S4)[1]
	}
	for(i in 1:length(diseases)){
		chars = apply(S4[subset,feature], 2, function(x) tapply(x, INDEX = as.factor(S4[subset, diseases[i]]), mean, na.rm= T))
		write.csv(t(chars), file = paste(names(diseases)[i],".csv", sep = ""))
	}
}

# partial correlation calculation
cor.partial<-function(data, variables, ...){
	P=cor(data[,variables], ...)
	P.invers=solve(P)
	P.diag = diag(P.invers)
	Z = -P.invers/sqrt(P.diag%*%t(P.diag))
	return(Z)
}

#correlation matrix to linkage pairs
cor2link = function(Zeta, threshold){
	ggm=NULL;
	for(i in 1:(dim(Zeta)[1]-1)){
		for(j in (i+1):dim(Zeta)[1]){
			if(abs(Zeta[i,j])>= threshold){
				tmp=c(colnames(Zeta)[i], colnames(Zeta)[j], 
						Zeta[i,j])
				ggm=rbind(ggm,tmp)
			}
		}
	}
	return(ggm)
}

#	Differential correxpression
#	The method was described in Cho, Sung, Jihun Kim, and Ju Kim. “Identifying Set-wise Differential Co-expression in Gene Expression Microarray Data.” BMC Bioinformatics 10, no. 1 (2009): 109.
#	Description of parameters:
#	cor1, cor2: two pearson correlation values
#	N1, N2: the number of samples used for the calculation of cor1 and cor2
diffcorr <- function(cor1, cor2, N1, N2){
	Z1=0.5*log((1+cor1)/(1-cor1))
	Z2=0.5*log((1+cor2)/(1-cor2))
	p=1-pnorm(abs((Z1-Z2)/sqrt(1/(N1-3)+1/(N2-3))))
	return(p)
}



