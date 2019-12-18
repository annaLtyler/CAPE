#This function pulls a threshold out of the data object
#based on a function supplied, ex. min, max, mean


get.threshold <- function(data.obj, fun){

	get.fun <- match.fun(fun)	
	
	all.threshold.locale <- grep("thresh", names(data.obj))
	all.fun.calls.locale <- grep("call", names(data.obj))
	threshold.calls.locale <- intersect(all.threshold.locale, all.fun.calls.locale)
	threshold.vals.locale <- setdiff(all.threshold.locale, threshold.calls.locale)
	threshold.vals <- apply(matrix(threshold.vals.locale, ncol = 1), 1, function(x) data.obj[[x]])

	final.threshold <- get.fun(threshold.vals)
	return(final.threshold)
	
	
}