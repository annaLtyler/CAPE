#This function samples the individuals in the genotype and cross objects
#with replacement for bootstrapping.
#perc.ind specifies what percentage of the original population size
#the sampled population size will be.


sample.for.bootstrap <- function(data.obj, geno.obj, perc.ind = 80){
	

	n.ind <- nrow(data.obj$pheno)
		
	n.fixed <- round((n.ind*perc.ind)/100)
	the.rest <- n.ind - n.fixed
	sampled.ind <- sample(1:n.ind, n.fixed, replace = FALSE)
	sampled.ind <- c(sampled.ind, sample(sampled.ind, the.rest, TRUE))
	
	data.obj$pheno <- data.obj$pheno[sampled.ind,]
	data.obj$geno.names$mouse <- data.obj$geno.names$mouse[sampled.ind]
	
	if(!is.null(data.obj$p.covar.table)){
		data.obj$p.covar.table <- data.obj$p.covar.table[sampled.ind,]
		}
	
	geno <- geno[sampled.ind,,]

	final.results <- list("data.obj" = data.obj, "geno.obj" = geno)
	return(final.results)
	
	
}