# Population characteristics
# 
# Author: tao.xu
###############################################################################


###data read-in
S4 = read.csv("data/K9512_Wang_Sattler_S4_tra300712.csv", sep = ";")
F4 = read.csv("data/K9512_Wang_Sattler_F4_tra300712.csv", sep = ";")
S4_F4_map = read.csv("data/K9512_Wang_S_S4_F4trans300712.csv", sep = ";")

####	additional data
add = read.csv("data/sas data/pv95_12_add.csv", colClass = "character")
add = apply(add, 2, as.numeric)
tmp = merge(S4, add, by.x = "ZZ_nr.x", by.y = "ZZ_nr")


#S4_valid_measure = 
#F4_valid_measure = 
load("F:/nicotine/F4_valid_measures.RData")

summary(S4[,1:31])
#missing values
missing_value = apply(S4[,1:31], 2, function(x) length(which (is.na(x))))
#samples with missing values
sample_availiable = apply(S4[,c(1:4, 6:31)], 1, function(x) length(which(is.na(x))))
which(sample_availiable !=0)
#

S4$waist2hip = S4$lttumf / S4$lthumf
F4$waist2hip = F4$uttumf / F4$uthumf

#S4 features and diseases
feature.cont = scan(what = character())
ltalteru
ll_choln
ll_hdln
ltsysmm
ll_ldln
ll_hgb
ll_hbav
ll_trin
lh_crp
lttumf
lthumf
ltdiamm
waist2hip
ltbmi
ltalkkon

feature.disc = scan(what = character())
ltcigreg
ltrauchp
lcsex
ltv_pdmt



diseases = scan(what = character())
ltjnc7
ltschl
ltmi
lc044f_1

names(diseases) = c("Hypertension","Stroke","Myocardio","Heart failure")

#population characteristics
for(i in 1:length(diseases)){
	chars = apply(S4[,feature.cont], 2, function(x) tapply(x, INDEX = as.factor(S4[, diseases[i]]), mean, na.rm= T))
	write.csv(t(chars), file = paste(names(diseases)[i],"S4 cross sectional.csv", sep = " "))
}


for(i in 1:4){
	print(diseases[i])
	print(tapply(S4$ltv_pdmt, INDEX = as.factor(S4[, diseases[i]]), function(x) length(which(x == 2))/ length(x) ))
}

################# F4 #################
F4$waist2hip = F4$uttumf/F4$uthumf

###additional features of F4:
#	1.utmsidfup: defintion of Metabolic Syndrom by International Diabetes Federation 2009, which be use as a indicator of high risk for CVD events
#	2.utmialt and utschalt: age at which first diagenosed for Myocardial infarction and stroke, respectively.
#	3.utsmkreg: parameter describing smoking behavior. differed from utcigreg that classified former smoker by their smoking habits (regular/occasional)
#	4.us_c07a: parameter describing the cardiac arrhythmias in the last 12 months
feature.cont = scan(what = character())
utalteru
ul_choln
ul_hdln
utsysmm
ul_ldln
ul_hgb
ul_hbav
ul_trin
uh_crp
uttumf
uthumf
utdiamm
waist2hip
utbmi
utalkkon

feature.disc = scan(what = character())
utcigreg
utrauchp
ucsex
utv_pdmt
utmsidfup

#ltv_pdmt, not in F4?? utv_pdmt is not applied in PV  

diseases = scan(what = character())
utjnc7
utschl
utmi
us_c04a

names(diseases) = c("Hypertension","Stroke","Myocardio","Heart failure")

#population characteristics
for(i in 1:length(diseases)){
	chars = apply(F4[,feature], 2, function(x) tapply(x, INDEX = as.factor(F4[, diseases[i]]), mean, na.rm= T))
	write.csv(t(chars), file = paste(names(diseases)[i],"F4 cross sectional.csv", sep = " "))
}


for(i in 1:4){
	print(diseases[i])
	print(tapply(F4$utmsidfup, INDEX = as.factor(F4[, diseases[i]]), function(x) length(which(x == 1)) ))
}

#####################	S4 F4 mapping	######################
index.prospect = apply(S4_F4_map, 1, function(x) !is.na(x[3])&!is.na(x[6]) ) 
S4_F4_map.prospect = S4_F4_map [index.prospect, ]
Cohort = data.frame(S4 = S4_F4_map.prospect[3], F4= S4_F4_map.prospect[6])
mode(Cohort)

Cohort = Cohort[which(as.character(Cohort[,1]) %in% rownames(S4)),] 

H = data.frame( s4.hf = S4[as.character(Cohort[,1]), "lc044f_1"], f4.hf = F4[as.character(Cohort[,2]),"us_c04a"])
S = data.frame ( S4 = S4 [as.character(Cohort[,1]),"ltjnc7"], F4 = F4[as.character(Cohort[,2]),"utjnc7"])
remove(H); remove(S)

######################	population characteristics of prospective data	#####################
#S4
source("D:/Users/tao.xu/Documents/GitHub/CVD/CVD/Toolkits.R")
for(i in 1:length(diseases.s4)){
	
	tmps4 = as.factor(S4$prev_apo)
	tmpf4 = as.factor(S4$inz_apo)
	
	chars = characteristics(S4[, feature.cont], factor = interaction(tmpf4, tmps4), d = 2, na.rm = T)
#	chars = apply(
#			S4[Cohort$zz_nr_s4,feature.cont], 2, 
#			function(x) {
#				m = tapply(x, INDEX = interaction(tmpf4, tmps4), mean, na.rm = T)
#				s = tapply(x, INDEX = interaction(tmpf4, tmps4), sd, na.rm = T)
#				m = round(m, digits = 1)
#				s = round(s, digits = 1)
#				return(paste(m, "(", s, ")", sep = ""))
#			} 
#	)
	write.csv(t(chars), file = "Stroke S4 prospective.csv")
}



for(i in 1:3){
	print(diseases[i])

	tmps4 = as.factor(S4$prev_mi)
	tmpf4 = as.factor(S4$inz_mi)
	#which(S4$lthyact==2)
	print(
			tapply(S4[, "lcsex"], 
					INDEX = interaction(tmpf4, tmps4), 
					function(x) length(which(x == 1))#/length(x)
			)
	)

}

apply(S4[which(S4$lthyact==1), feature.cont], 2, function(x) wilcox.test(x~S4$inz_apo[which(S4$lthyact==1)]))

#F4

for(i in 1:length(diseases)){
	chars = apply(F4[as.character(Cohort[,2]),feature.cont], 2, function(x) tapply(x, INDEX = as.factor(F4[as.character(Cohort[,2]), diseases[i]]), mean, na.rm= T))
	write.csv(t(chars), file = paste(names(diseases)[i],"F4 prospective.csv", sep = " "))
}

for(i in 1:4){
	print(diseases[i])
	print(tapply(F4[as.character(Cohort[,2]),"utmsidfup"], INDEX = as.factor(F4[as.character(Cohort[,2]), diseases[i]]), function(x) length(which(x == 1)) ))
}



hist(S4 [as.character(Cohort[,1]),"ltsysmm"],breaks = 13, col = "blue", ylim = c(0,300), xlim = c(70, 200), main = "Histogram of systolic blood pressure in S4(blue) and F4 (red)", ylab = "Frequency", xlab = "Systolic blood pressure")
hist(F4 [as.character(Cohort[,2]),"utsysmm"], breaks = 13, col = "red", add=T, xlim = c(70, 200))

hist(S4 [as.character(Cohort[,1]),"ltdiamm"],breaks = 12, col = "blue", ylim = c(0,250), xlim = c(30, 150), main = "Histogram of diatolic blood pressure in S4(blue) and F4 (red)", ylab = "Frequency", xlab = "Diatolic blood pressure")
hist(F4 [as.character(Cohort[,2]),"utdiamm"], breaks = 12, col = "red", add=T, xlim = c(30, 150))


# age check
tmp = data.frame ( S4 = S4 [as.character(Cohort[,1]),"ltalteru"], F4 = F4[as.character(Cohort[,2]),"utalteru"])



####################################################
tmp = read.csv("F:/nicotine/Diabetes and Smoking/S4_addition.csv")
tmp = tmp[,-5]
tmp = merge(tmp, S4)