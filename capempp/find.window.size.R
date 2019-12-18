#This function finds a window size for binning a genotype
#object into chunks for correlation testing
#it operates on one chromosome at a time

find.window.size <- function(data.obj, geno.obj, chr = 1, min.cor = 0.2){

	#================================================================
	#internal functions
	#================================================================
	get.cor <- function(marker.idx){
		marker1.geno <- as.vector(geno.obj[,,marker.idx[1]])
		marker2.geno <- as.vector(geno.obj[,,marker.idx[2]])
		return(cor(marker1.geno, marker2.geno))
		}
	#================================================================

	
	all.ch <- data.obj$chromosome
	
	ch.locale <- which(all.ch == chr)

	start.marker <- 1
	next.marker <- 2
	marker.cor <- get.cor(c(start.marker, next.marker))
	all.num.to.min <- NULL
	
	while(next.marker < length(ch.locale)){
		while(marker.cor > min.cor){
			next.marker <- next.marker + 1
			marker.cor <- get.cor(c(start.marker, next.marker))		
			}		
		all.num.to.min <- c(all.num.to.min, next.marker)
		start.marker <- next.marker
		next.marker <- start.marker + 1
		marker.cor <- get.cor(c(start.marker, next.marker))
		}
	
	pos.pairs <- consec.pairs(all.num.to.min)
	window.sizes <- apply(pos.pairs, 1, function(x) x[2] - x[1])

	return(max(window.sizes))
	
	
	
	
}