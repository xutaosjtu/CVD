# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


Name = names(diseases.s4)
Models = c("Model 1", "Model 2", "Model 3")
dir("Cross sectional\\")

for (n in Name){
	for(j in 1:3){
		file1 = paste("Cross sectional/", "Model " , j , "/",  n, " categorical model_", j, " S4.csv",  collapse = "", sep ="")
		file2 = paste("Cross sectional/", "Model " , j , "/",  n, " categorical model_", j, " F4.csv",  collapse = "", sep ="")
		rst.s4 = read.csv(file1);rst.f4 = read.csv(file2)
		rst = merge(rst.s4, rst.f4, by.x = "X", by.y = "X", all = TRUE)
		file3 = paste("Cross sectional/", "Model " , j , "/",  n, " categorical model_", j, ".csv",  collapse = "", sep ="")
		write.csv(rst, file3, row.names = F)
	}
}


rst.s4 = read.csv("F:/Cardiovascular disease/Cross sectional/Model 4/CVD combined_Model SCORE_S4.csv")
rst.f4 = read.csv("F:/Cardiovascular disease/Cross sectional/Model 4/CVD combined_Model SCORE_F4.csv")
rst = rst = merge(rst.s4, rst.f4, by.x = "X", by.y = "X", all = TRUE)
write.csv(rst, file = "Cross sectional/CVD combined_Model SCORE.csv")