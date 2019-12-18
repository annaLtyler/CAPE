#This function finds diagonal indices of a matrix
#starting from the x = y diagonal. This is for 
#identifying pairs of adjacent elements for calculating
#correlations between
#depth indicates how many diagonals you want to include
#in the indices

pair.idx <- function(idx, depth = 1){
		
	mat1 <- matrix(idx, length(idx), length(idx), byrow = TRUE)
	mat2 <- matrix(idx, length(idx), length(idx), byrow = FALSE)
	
	diff.mat <- mat1 - mat2
	diff.mat[which(diff.mat < 0)] <- NA
	
	all.idx	<- which(diff.mat <= depth, arr.ind = TRUE)
	
	return(all.idx)
	} #end function