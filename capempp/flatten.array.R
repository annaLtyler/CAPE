#This function takes in a 3D array of values, effects or t.stats, etc.
#and returns a 2D matrix in wich each entry contains a single value 
#representing all corresponding entries in the 3rd dimension. This 
#number is determined by the function the user uses. It could be min, 
#max, mean, etc.

flatten.array <- function(arrayX, margin1, margin2, slice.fun){
	
	dims <- 1:length(dim(arrayX))
	array.dim <- unlist(dim(arrayX))
	final.row <- array.dim[margin1]
	final.col <- array.dim[margin2]
	
	slice.dim <- array.dim[-margin1]
	new.margin2 <- which(slice.dim == final.col)
	
	apply.fun <- match.fun(slice.fun)

	flattened.mat <- apply(arrayX, margin1, function(slice) apply(slice, new.margin2, apply.fun))

	return(flattened.mat)


}