# TODO: Add comment
# 
# Author: tao.xu
###############################################################################

##data
tmp = data.frame(
		time = S4$mi_time, event = S4$inz_mi,  # time and events
		S4$ltalteru, S4$ltbmi, as.factor(S4$lcsex), ##model1
		(S4$lp_diab_who06==4|S4$lp_diab_who06==5), ##model 2
		log(S4$ltsysmm),log(S4$ll_hdln), log(S4$ll_choln), as.factor(S4$ltcigreg),##model 3
		##log(S4$total2HDL), ##model 4
		##log(S4$ltdiamm), ##model 5
		log(as.matrix(S4[, metabo.asso]))
)
colnames(tmp)[3:10] = clinical#, "total2HDL"
na.index = unique(unlist(apply(tmp, 2, function(x) which(is.na(x)))))
subset = setdiff(which(S4$prev_mi == 0 & !is.na(S4$inz_mi)), na.index)

require(CoxBoost)
# find optimal penalty
penalty.optimal = optimCoxBoostPenalty(
		time = S4$mi_time[subset], 
		status = S4$inz_mi[subset],
		x = as.matrix(tmp[subset, 12:31]),
		minstepno = 1, 
		maxstepno = 100)
# find optimal boosting steps
booststeps.optimal = cv.CoxBoost(
		time = S4$mi_time[subset],
		status = S4$inz_mi[subset], 
		x = log(as.matrix(S4[subset, c(metabolites.asso)])),  
		type="verweij")
#variable selection
MI.metabocox = CoxBoost(
		time = S4$mi_time[subset ], 
		status = S4$inz_mi[subset ], 
		#unpen.index = c(3:9),
		x =  as.matrix(tmp[subset, c(metabo.asso, clinical)]),
		stepno = 89, 
		penalty = 603)

metabo.selected = scan(what = character())
C14_1
C16_1
Arg
Trp
lysoPC_a_C17_0
lysoPC_a_C18_2
PC_aa_C28_1
PC_aa_C32_2
PC_ae_C36_2
PC_ae_C38_0
PC_ae_C40_1

###evaluation by ROC
require(ROCR)
pred = prediction( 1-predict(MI.metabocox, type = "risk", times = 3683), tmp$event[subset])
plot(performance(pred, "tpr", "fpr"), col = 5)
performance(pred, "auc")
 
#by boost in gbm which uses gradient descent method
require(gbm) 
MI.gbmboost = gbm(formula = y ~ .,
		distribution = "coxph",
		data = as.data.frame(tmp[subset,])
)
data = data.frame( S4$mi_time, S4$inz_mi, tmp)
