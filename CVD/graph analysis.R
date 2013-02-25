# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


library(igraph)
g=read.graph("KEGG Reactions.csv",format="ncol",directed=T)

bcc=betweenness(g)
names(bcc)=V(g)$name
degree=degree(g)
names(degree)=V(g)$name
cc=closeness(g)
names(cc)=V(g)$name

metabolites=data.frame(
		metabolites
		,bcc=bcc[metabolites]
		,degree=degree[metabolites]
		,cc=cc[metabolites]
		)

metabolites=metabolites[complete.cases(metabolites),]
		
plot(density(log(metabolites$bcc),adjust=0.7),col="red")
lines(density(log(bcc),adjust=0.7))
		
plot(density(log(degree),adjust=0.9))
lines(density(log(metabolites$degree),adjust=0.9),col="red")

####################	Find shortest path	##########################
require(igraph)
score = 0.9

all=read.csv('F:/Database/meta_enzyme_string_unique.v9.0',head=F, sep = "\t")
all=all[all[,3] >= score,]

graph=read.graph('F:/Database/meta_enzyme_string_unique.v9.0',format='ncol')
del.edge=E(graph)[E(graph)$weight < score]
graph=delete.edges(graph, del.edge)

metabolites = as.character(unique(all[all[,3]==1,1]))
#meta.id = which(V(graph)$name %in% metabolites)-1
enzymes = as.character(unique(all[all[,3]==1,2]))
enzymes = as.character(unique(all[which(all[,1] %in% From),2]))

From =scan(what = character())
Arg-PTC
Trp-PTC
PC-aa-C32:2
lysoPC-a-C17:0

SM-C24:0
PC-ae-C38:6
PC-aa-C32:1
SM-OH-C22:2
His-PTC

To=c(unique(as.matrix(read.csv('F:/Cardiovascular disease/SNPs/CAD gene GWAS.csv',head=F))),
		unique(as.matrix(read.csv('F:/Cardiovascular disease/SNPs/CAD gene FDR.csv', header = F)))
)

tmp=unique(c(as.character(all[,1]),as.character(all[,2])))
setdiff(To,tmp)
To=intersect(To,tmp)

fid = which(V(graph)$name %in% From)
tid = which(V(graph)$name %in% To)

path2link = function(path){
	link = NULL;
	for(i in 1:(length(path)-1)){
		link = rbind(link, c(path[i], path[i+1]))
	}
	return(link)
}

network = NULL;
for (i in 1:length(fid)){
	path = get.shortest.paths(graph, from=fid[i], to=tid,
			weights = E(graph)$weight, output = "vpath"
	)
	
	path = sapply(path, function(x) V(graph)$name[x])
	if(i ==2){
		path = path[which(sapply(path, length)<=4)]	
	}
	else path = path[which(sapply(path, length)<=4)]
	link = lapply(path, path2link)
	for(j in 1:length(link)){
		network = rbind(network, link[[j]])
	}
}


nodes = unique(c(enzymes, To, metabolites))
rst = data.frame(nodes, 
		'type' = rep("NA", length(nodes)),
		'GWAS' = rep("NA", length(nodes))
)

rst = as.matrix(rst)
rst[which(rst[,1] %in% enzymes), 2] = "enzymes"
rst[which(rst[,1] %in% metabolites), 2] = "metabolites"
rst[which(rst[,1] %in% To), 3] = "GWAS"


meta_enzyme = all[which(all[,1] %in% From),]

write.csv(meta_enzyme, file = "F:/Cardiovascular disease/SNPs/metabo enzyme paris.csv")

To[which(To == "AYTL2")] = "LPCAT1"
To[which(To == "AYTL1")] = "LPCAT2"
To[which(To == "TMEM23")] = "SGMS1"
To[which(To == "INDO")] = "IDO1"
To[which(To == "NOS2A")] = "NOS2"
To[which(To == "MBOAT5")] = "LPCAT3"
To[which(To == "APO5A")] = "APOA5"
