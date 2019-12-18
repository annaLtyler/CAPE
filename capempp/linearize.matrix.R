#This script takes in a matrix with row and columnames and returns
#a linearized version of the matrix with the appropriate names for
#each entry. The user can specify the upper triangle or the lower
#triangle of the matrix and whether they want the diagonal entries
#or not.

linearize.matrix <- function(mat, upper.triangle = TRUE, diagonal = TRUE){
	
	mat.dim <- dim(mat)
	if(mat.dim[1] != mat.dim[2]){
		stop("The matrix must be square")
		}
	
	row.mat <- matrix(1:mat.dim[1], mat.dim[1], mat.dim[1], byrow = FALSE)
	col.mat <- matrix(1:mat.dim[1], mat.dim[1], mat.dim[1], byrow = TRUE)
	
	if(upper.triangle){
		rowname.indices <- row.mat[upper.tri(row.mat, diag = diagonal)]
		colname.indices <- col.mat[upper.tri(col.mat, diag = diagonal)]
		
		all.names <- NULL
		for(i in 1:length(rowname.indices)){
			all.names <- c(all.names, paste(rownames(mat)[rowname.indices[i]], colnames(mat)[colname.indices[i]], sep = "_"))
		}
	
		result <- mat[upper.tri(mat, diag = diagonal)]
		names(result) <- all.names	
		return(result)

	}else{
		
		rowname.indices <- row.mat[lower.tri(row.mat, diag = diagonal)]
		colname.indices <- col.mat[lower.tri(col.mat, diag = diagonal)]
		
		all.names <- NULL
		for(i in 1:length(rowname.indices)){
			all.names <- c(all.names, paste(rownames(mat)[rowname.indices[i]], colnames(mat)[colname.indices[i]], sep = "_"))
		}
	
		result <- mat[lower.tri(mat, diag = diagonal)]
		names(result) <- all.names	
		return(result)	
		}
	
}