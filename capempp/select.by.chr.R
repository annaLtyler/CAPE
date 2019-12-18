#This script subsets a cross by the given
#chromosome number.
#The cross object is returned with only 
#genotype data for the selected chromosomes
#just modify the data.obj, we will only scan
#based on the names in the data.obj

select.by.chr <- function(data.obj, chr){
		
	chr.list <- data.obj$chromosome
	chr.locale <- which(chr.list %in% chr)
	
	data.obj$chromosome <- data.obj$chromosome[chr.locale]
	data.obj$marker.location <- data.obj$marker.location[chr.locale]
	data.obj$marker.num <- data.obj$marker.num[chr.locale]
	data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][chr.locale]
	
	return(data.obj)
	
	}