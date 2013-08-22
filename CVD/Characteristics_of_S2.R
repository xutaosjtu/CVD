S2$total2HDL = S2$cl_chola/S2$cl_hdla
S2$Arg.Trp = S2$Arg/S2$Trp

data.S2 = subset(S2, subcoho==1 | inz_mi02==1)
data.S2 = subset(data.S2, prev_mi==0)

feature.cont = scan(what = character())
ctalteru
ctbmi
ctalkkon
cl_chola
cl_hdla
cl_ldla
total2HDL
ctdiamm
ctsysmm
cl_crp
cttumf
cthumf
waist2hip

chars = characteristics(data.S2[, feature.cont], factor = as.factor(data.S2$inz_mi02), d = 2, na.rm = T)
sapply(data.S2[,feature.cont], function(x) wilcox.test(x~data.S2$inz_mi)$p.value)

tapply(data.S2$ccsex, INDEX = data.S2$inz_mi02, function(x) table(x)/length(x))
fisher.test(table(data.S2$ccsex, data.S2$inz_mi02))
