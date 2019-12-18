#This script subsets a cross by the given
#chromosome number.
#The cross object is returned with only 
#genotype data for the selected chromosomes

subset.by.chr <- function(data.obj, chr){
	
	geno <- data.obj$geno
	
	#Get the dimension names to minimize confusion	
	mouse.dim <- which(names(dimnames(geno)) == "mouse")
	locus.dim <- which(names(dimnames(geno)) == "locus")
	allele.dim <- which(names(dimnames(geno)) == "allele")
	
	chr.list <- data.obj$chromosome
	chr.locale <- which(chr.list %in% chr)
			
	if(length(dim(geno)) == 3){
		sub.geno <- geno[,,chr.locale]
		}

	names(dimnames(sub.geno)) <- names(dimnames(geno))
	data.obj$geno <- sub.geno
	data.obj$chromosome <- chr.list[chr.locale]
	data.obj$marker.location <- data.obj$marker.location[chr.locale]

	return(data.obj)
	
	}