setwd("../../../Dropbox/Cardiovascular disease/")
data<-read.csv("data/K9512_Wang_Sattler_S2.csv")

tmp = grep("X", colnames(data), fixed=T)
data = data[, -tmp]
measures = colnames(data)[19:249]
tmp = sapply(measures, function(x) tapply(data[,x], INDEX = data$Sample.Bar.Code, mean, na.rm=T))

###################################################
# Overall quality of the measurement
# 1. Overall CV and within plate CV
# 2. Between plate differences of refsamples
# 3. Measurements above LOD
###################################################
index.ref = sapply(data$Sample.Identification, function(x) grep("Ref",x,fixed=T) )
index.ref = sapply(index.ref, function(x) length(x)!=0)

##over all CV and within plate CV 
rst=NULL
for (i in measures){
  within = tapply(data[index.ref,i], 
                  INDEX = data$Plate.Bar.Code[index.ref], 
                  function(y) sd(y,na.rm=T)/mean(y,na.rm=T)
  )
  overall = sd(data[index.ref,i],na.rm=T)/mean(data[index.ref,i],na.rm=T)
  rsti = c(within,overall)
  rst = cbind(rst,rsti)
}
colnames(rst) = measures
write.csv(t(rst), file = "CV of each metabolites across all paltes.csv")

## between plate differences of refsamples 
pdf("reference samples in different plates.pdf", width=15, height=15)
par(mfrow = c(5,5))
for(i in measures){
  data$m = data[,i]
  #plotmeans(m ~ Plate.Bar.Code, data, subset=index.ref)
  boxplot(m ~ Plate.Bar.Code,data,subset=index.ref, main = i)
}
dev.off()

## Measurements above LOD
aboveLOD = function(data, measures){
  # measurement above LOD in one plate
  # input: meatbolite measurements of samples and zero samples in one plate
  # output: Logical indicator of above LOD
  index.zero = which(data$Sample.Type=="Zero Sample")
  index.sample = grep("30", data$Sample.Identification, fixed=T)
  #index.nurse = sapply(index.nurse, function(x) length(x)!=0)
  #print(head(data))
  rst.overLOD=sapply(data[,measures], 
                    function(x) {
                      medi=median(x[index.zero], na.rm=T)
                      #print(medi)
                      if(!is.na(medi)) return(x>3*medi)
                      else return(rep(NA,length(x)))
                    }
                    
  )
  return(rst.overLOD)
}

matrixLOD = data[,measures] # indicator matrix of measurements above LOD
for(i in names(table(data$Plate.Bar.Code))){
  plate = which(data$Plate.Bar.Code == i)
  matrixLOD[plate,]=aboveLOD(data[plate,], measures)
}
matrixLOD = data.frame(matrixLOD, 
                       Sample.Identification=data$Sample.Identification)
matrixLOD = merge(matrixLOD,samples, by.x = "Sample.Identification", by.y = "Proben_ID")

rst=NULL # measurements above LOD in each plate and all plates
for(i in names(table(data$Plate.Bar.Code))){
  subset = which(data$Plate.Bar.Code == i)
  tmp = aboveLOD(data[subset,])
  rst = cbind(rst,apply(tmp,2,sum, na.rm = T))
}
colnames(rst) = names(table(data$Plate.Bar.Code))
rst = data.frame(rst, overall = apply(rst,1,sum))
rst$overall=rst$overall/nrow(samples)
write.csv(rst, file = "overLOD_nurse.csv")

## data Normalization
data.merged = merge(data,samples,by.x="Sample.Identification", by.y="Proben_ID")

normalize<-function(data, measures){
  tmp = apply(data[,measures], 2, function(x) 1000*x/data$Creatinine)
  data[,measures] = tmp
  return(data)
}

data.merged = data.merged[which(matrixLOD$Creatinine==1), ]# exclude abnormal creatinine data
matrixLOD = matrixLOD[which(matrixLOD$Creatinine==1),]# exclude abnormal creatinine data
data.merged = normalize(data.merged,measures)

for(i in measures){
  data.merged[which(matrixLOD[,i]==0),i]=NA
}

## outliers
index=apply(data.merged[,measures],2, function(x) which(abs(x)>mean(x,na.rm=T)+5*sd(x,na.rm=T)|abs(x)<mean(x,na.rm=T)-5*sd(x,na.rm=T))) 
for(i in measures){
  data.merged[index[[i]],i]=NA
}
