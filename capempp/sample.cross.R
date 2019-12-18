#This function samples the individuals in the genotype and cross objects
#with replacement for bootstrapping.
#if n.ind is NULL, the sample will be the same number as individuals


sample.cross <- function(data.obj, geno.obj, n.ind = NULL, with.replacement = TRUE){
	
	if(is.null(n.ind)){
		n.ind <- nrow(data.obj$pheno)
		}
		
	sampled.ind <- sample(1:n.ind, n.ind, replace = with.replacement)
	data.obj$pheno <- cross$pheno[sampled.ind,]
	data.obj$geno.names$mouse <- data.obj$geno.names$mouse[sampled.ind]
	
	if(!is.null(data.obj$p.covar.table)){
		data.obj$p.covar.table <- data.obj$p.covar.table[sampled.ind,]
		}
	
	geno <- geno[sampled.ind,,]

	final.results <- list("data.obj" = data.obj, "geno.obj" = geno)
	return(final.results)
	
	
}