#This function quickly grabs individual 
#values from a matrix based on a vector
#of row indices and a vector of column
#indices

mat.vals <- function(matX, row.ind, col.ind){		
	not.na.locale <- intersect(which(!is.na(row.ind)), which(!is.na(col.ind)))
	
	#either of these methods works, but the first, is faster, I think
	vals <- matX[row.ind[not.na.locale]+nrow(matX)*(col.ind[not.na.locale] - 1)]
	# vals <- matX[cbind(row.ind[not.na.locale], col.ind[not.na.locale])]
	final.mat <- cbind(row.ind[not.na.locale], col.ind[not.na.locale], vals)	
	return(final.mat)
	
	}