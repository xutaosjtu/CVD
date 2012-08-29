# TODO: Add comment
# 
# Author: tao.xu
###############################################################################


Name = names(diseases)
Models = c("Model 1", "Model 2", "Model 3")
dir("Cross sectional\\")

for(j in 1:3){
	for (n in Name){
		file1 = paste("Cross sectional/", "Model " , j , "/",  n, " model_", j, " S4.csv",  collapse = "", sep ="")
		file2 = paste("Cross sectional/", "Model " , j , "/",  n, " model_", j, " F4.csv",  collapse = "", sep ="")
		rst.s4 = read.csv(file1);rst.f4 = read.csv(file2)
		rst = merge(rst.s4, rst.f4, by.x = "X", by.y = "X", all = TRUE)
		file3 = paste("Cross sectional/", "Model " , j , "/",  n, " model_", j, ".csv",  collapse = "", sep ="")
		write.csv(rst, file3, row.names = F)
	}
}
