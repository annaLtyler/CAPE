#This function subsets a cape data.obj for testing purposes

subsample.cape <- function(data.obj, n.markers = 1000){
	sampled.snps <- sort(sample(1:length(data.obj$geno.names[[3]]), 1000))
	data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][sampled.snps]
	data.obj$marker.num <- data.obj$marker.num[sampled.snps]
	data.obj$chromosome <- data.obj$chromosome[sampled.snps]
	data.obj$marker.location <- data.obj$marker.location[sampled.snps]
	return(data.obj)
	}