#This function finds the t statistic threshold based on 
#permutations and using the extreme value distribution

p.val.thresh <- function(perm.max, alpha){
	
	require("evd")

	get.s <- function(evd.result){
		s <- qgev(1-alpha,loc=evd.result$estimate[1], scale=evd.result$estimate[2], shape=evd.result$estimate[3], lower.tail = TRUE)
		return(s)
		}


	#apply the extreme value distribution to the results
	evd <- apply(perm.max, 2, function(x) fgev(x, std.err = FALSE))

	
	s <- as.vector(sapply(evd, get.s))

	#calculate one threshold over all phenotypes
	threshold <- mean(s)

	return(threshold)	
	
}