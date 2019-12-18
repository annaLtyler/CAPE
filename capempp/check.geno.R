	#check to see if the genotype matrix is one of the covariates
	
		check.geno <- function(genotype, covar.table){

			if(is.null(covar.table)){
				return(NULL)
				}
			
			covar.which <- NULL
			for(i in 1:ncol(covar.table)){
				
				if(identical(as.vector(covar.table[,i]), as.vector(genotype[,1]))){
					covar.which <- c(covar.which, i)
					}
				}
				
			return(covar.which)
			}
