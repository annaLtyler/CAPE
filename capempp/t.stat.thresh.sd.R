#This function takes in t statistics from wrapper scanone and
#calculates a threshold based on the mean and standard deviation


t.stat.thresh.sd <- function(t.stats, num.sd = 2){
	
	
	mean.t <- mean(t.stats, na.rm = TRUE)
	sd.t <- sd(t.stats, na.rm = TRUE)

	threshold <- mean.t + (num.sd*sd.t)
	return(threshold)
	
	
}