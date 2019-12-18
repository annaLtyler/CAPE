#This function recalculates a genomewide threshold for
#either the pairscan threshold or the covaraite threshold
#based on statistics that can be calculated from the 
#t statistic matrix

#it also adjusts the covariate flags accordingly


recalculate.threshold <- function(data.obj, pairscan.thresh.fun = NULL, pairscan.thresh.param = NULL, covar.thresh.fun = NULL, covar.thresh.param = NULL){
	
	t.vals <- data.obj$oneDscan.t.stats


	if(!is.null(pairscan.thresh.fun)){
		fun <- match.fun(pairscan.thresh.fun)
		new.threshold <- fun(t.vals, as.numeric(pairscan.thresh.param))
	
		data.obj$pairscan.thresh <- new.threshold
		data.obj$pairscan.thresh.call <- paste(pairscan.thresh.fun, "param =", pairscan.thresh.param)
		}
		

	if(!is.null(covar.thresh.fun)){
		fun <- match.fun(covar.thresh.fun)
		new.threshold <- fun(t.vals, as.numeric(covar.thresh.param))
	
		data.obj$covar.thresh <- new.threshold
		data.obj$covar.thresh.call <- paste(covar.thresh.fun, "param =", covar.thresh.param)

		#adjust the covariate flags. Each allele is a covariate if 
		#it's t statistic exceeds the covar threshold
		
		new.covar <- get.significant.loci(data.obj, new.threshold)
		covar.flags <- data.obj$covar.flags
		new.covar.flags <- array(0, dim = dim(covar.flags))
		dimnames(new.covar.flags) <- dimnames(covar.flags)
		#for each of the new covariates, find their location in the 
		#flags array and set the value to 1
		for(cv in 1:length(new.covar[,1])){
			locus <- new.covar[cv,"locus"]
			locus.locale <- which(dimnames(new.covar.flags)[[1]] == locus)
			pheno <- new.covar[cv,"phenotype"]
			pheno.locale <- which(dimnames(new.covar.flags)[[2]] == pheno)
			allele <- new.covar[cv,"allele"]
			allele.locale <- which(dimnames(new.covar.flags)[[3]] == allele)
			new.covar.flags[locus.locale, pheno.locale, allele.locale] <- 1
			}
		data.obj$covar.flags <- new.covar.flags


		} #end case for adjusting the covar.thresh

		

	return(data.obj)
	
}