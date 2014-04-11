## By a matched case control study
data = S4[which(S4$prev_mi==0),c("ltalteru", "lcsex", "ltbmi","my.diab","ltsysmm","my.cigreg","my.alkkon","ll_chola", "ll_hdla", "ll_ldla", "lh_crp", "inz_mi", "mi_time", S4_valid_measures)]
data = na.omit(data)
data$my.cigreg = as.numeric(data$my.cigreg)
data$my.diab = as.numeric(data$my.diab)
data$Arg_Trp = data$Arg/data$Trp

#Y = 
require(Matching)
tr = data$inz_mi
x = data[,c("ltalteru", "lcsex")]
set.seed(27)
rst.match = Match(Tr = tr, X =x, exact = T, M = 2, ties = F)

#rst.balance = MatchBalance(inz_mi ~ scale(ltalteru) + as.factor(lcsex) + scale(ltbmi)## model 1
#+ my.diab  ##model 2
#+ scale(ltsysmm) + my.cigreg + my.alkkon  + scale(ll_chola) + scale(ll_ldla) ##model 3+ total2HDL
#+ scale(lh_crp)  ##model 4
#, data, match = rst.match, nboots = 10)

data = data[c(rst.match$index.control, unique(rst.match$index.treated)), ]
#data$ID = rep(1:table(data$inz_mi)[1], 2)
data$ID = c(rep(1:table(data$inz_mi)[2], each = 2),1:table(data$inz_mi)[2])
require(survival)
rst = NULL
for(m in c(S4_valid_measures, "Arg_Trp")){
  #data$m = log(data[,"Arg"])-log(data$Trp)
  data$m = scale(log(data[,m]))
  model = clogit(inz_mi ~ m 
                #+ scale(ltbmi)## model 1
                #+ scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
                #+ as.factor(my.diab) + as.factor(my.cigreg) + as.factor(my.alkkon)
                #+ scale(log(lh_crp))  ##model 4
                + strata(ID)
                , data)
  rst = rbind(rst, summary(model)$coef[1,])
}
rownames(rst) = c(S4_valid_measures, "Arg_Trp")

set.seed(20)
rst = NULL
for(j in 1:40){
  #mstat = NULL;
  for(i in 1:100){
    tmp = analysis.matched(data, reduce = j, mnum=4, metabo = "Arg_Trp")
    rst= rbind(rst, tmp)
  }
  #z.est = c(mean(abs(rst[,4])), se = sd(abs(rst[,4]))/sqrt(40))
}
rst = data.frame(rst, index = rep(1:40, each = 100))

require(gplots)
plotmeans(z~index, data = rst, n.label=F,main = "Arginine tryptophan ratio")
#qnorm(0.025)
abline(h = qnorm(0.975), col = "red")
#1-pnorm(mean(abs(rst[,4])))

plot(-log10(rst[,5]))
abline(h = -log10(0.05))

#rst=boot(data, analysis.matched,R =100, mnum =4, metabo ="Arg")
#
#plot(-log10(rst$t[,5]))
#abline(h = -log10(0.05))
#
#sum(rst$t[,4]>rst$t0[4])/(1+rst$R)

analysis.matched <- function(data, reduce, mnum, metabo)
{
  require(survival)
  require(Matching)
  #print(indices)
  test = sample(which(data$inz_mi==1), reduce)
  train = data[-test,]
  tr = train$inz_mi
  x = train[,c("ltalteru", "lcsex", "my.diab","my.cigreg","my.alkkon")]
  rst.match = Match(Tr = tr, X =x, exact = T, M = mnum, ties = F)
  train = train[c(rst.match$index.control, unique(rst.match$index.treated)), ]
  train$ID = c(rep(1:table(train$inz_mi)[2], each = mnum),1:table(train$inz_mi)[2])
  train$m = log(train[,metabo])
  #train$m = log(train[,m])
  model = clogit(inz_mi ~ m 
                 + scale(ltbmi)## model 1
                 + scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla)##model 3
                 + scale(lh_crp)  ##model 4
                 + strata(ID)
                 ,train)
  return(summary(model)$coef[1,])
}


#data$ID = rep(1:table(data$inz_mi)[1], 2)

require(survival)
rst = NULL
for(m in S4_valid_measures){
  train$m = log(train$Arg)-log(train$Trp)
  #train$m = log(train[,m])
  model = clogit(inz_mi ~ m 
                 + scale(ltbmi)## model 1
                 + scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
                 + scale(lh_crp)  ##model 4
                 + strata(ID)
                 ,train)
  rst = rbind(rst, summary(model)$coef[1,])
}
rownames(rst) = S4_valid_measures

candidate = rownames(rst)[which(rst[,5]<0.05)]
train$Arg_Trp = train$Arg/train$Trp
model = clogit(inz_mi ~ Ala + Arg_Trp +Gly+Met+Phe+Ser+Trp+lysoPC_a_C16_0+lysoPC_a_C17_0+lysoPC_a_C18_0+lysoPC_a_C18_2+PC_aa_C28_1+PC_ae_C36_4+PC_ae_C36_5+PC_ae_C38_5
               + scale(ltbmi)## model 1
               + scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
               + scale(lh_crp)  ##model 4
               + strata(ID)
               ,train)

stepAIC(model)

model = clogit(inz_mi ~ Ala + Arg_Trp +Gly+lysoPC_a_C17_0+PC_ae_C36_5
               + scale(ltbmi)## model 1
               + scale(ltsysmm) + scale(ll_chola) + scale(ll_hdla) ##model 3+ total2HDL
               + scale(lh_crp)  ##model 4
               + strata(ID)
               ,train)

testset = data[c(test,which(data$inz_mi==0)),]
tr = testset$inz_mi
x = testset[,c("ltalteru", "lcsex", "my.diab","my.cigreg","my.alkkon")]
rst.match = Match(Tr = tr, X =x, exact = T, M = 1, ties = F)
testset = testset[c(rst.match$index.control, unique(rst.match$index.treated)), ]
testset$ID = c(rep(1:table(testset$inz_mi)[2], each = 1),1:table(testset$inz_mi)[2])

pred = predict(model, testset, type = "risk")
require(pROC)
plot(roc( testset$inz_mi, pred))