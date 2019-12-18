#This function gets pair indices for looking at 
#marker correlations on a chromosome
#It gets the indices of the pairwise elements
#along the diagonal of a square matrix up to the 
#depth requested.
#num is the number of elements you want to get
#pairs for, and depth is the number of rows out
#from the diagonal that you want to calculate

pair.idx <- function(num, depth, self.pairs = FALSE, directed = FALSE){
	
	mat1 <- matrix(1:num, nrow = num, ncol = num, byrow = TRUE)
	mat2 <- matrix(1:num, nrow = num, ncol = num, byrow = FALSE)
	mat.diff <- mat1 - mat2
	
	if(!self.pairs){
		mat.diff[which(mat.diff == 0)] <- NA
		}
	if(!directed){
		mat.diff[which(mat.diff < 0)] <- NA
		}
	
	the.pairs <- which(abs(mat.diff) <= depth, arr.ind = TRUE)
	
	# test.mat <- matrix(0, num, num)
	# test.mat[the.pairs] <- 1
	return(the.pairs)	
	
}