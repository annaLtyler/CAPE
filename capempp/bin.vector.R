#This function bins a continuously valued vector
#based on user-defined bins.
#it is useful for binning continuously valued genotypes
#each value in the matrix gets shifted to the nearest
#provided in the argument bins

bin.vector <- function(vectorX, bins = seq(0,1,0.5)){
	
	dist.mat <- apply(matrix(bins, ncol = 1), 1, function(x) x - vectorX)	
	binned.vector <- apply(dist.mat, 1, function(x) bins[which(abs(x) == min(abs(x)))[1]])
	return(binned.vector)	
	
}