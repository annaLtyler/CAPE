#This script allows users to change covariate assignments
#manually
# markers <- rownames(covar.flags)[5:10]
# pheno <- "ET1"
	


set.covar <- function(data.obj, singlescan.obj = NULL, markers = NULL, covar.thresh = NULL, pheno = NULL, allele = NULL, is.covar = TRUE, scan.what = c("eigentraits", "raw.trait")){
	
	check.underscore(data.obj)
	# check.bad.markers(data.obj, geno.obj)
	
	use.eigentraits <- length(c(grep("eigen", scan.what), grep("ET", scan.what), grep("et", scan.what)))
	if(use.eigentraits){
		cat("Setting covariates for eigentraits...\n")
		pheno.mat <- data.obj$ET
		}else{
			pheno.mat <- data.obj$pheno
			cat("Setting covariates for raw traits...\n")
			}

	
		t.stats <- singlescan.obj$singlescan.t.stats
	
		if(is.null(covar.thresh)){
			covar.thresh <- data.obj$covar.thresh
			}else{
				data.obj$covar.thresh <- covar.thresh
				}
					
	geno.names <- data.obj$geno.names
	
	covar.flags <- array(0, dim = c(length(geno.names[[3]]), dim(pheno.mat)[2], length(geno.names[[2]])))
	dimnames(covar.flags) <- list(geno.names[[3]], colnames(pheno.mat), geno.names[[2]])
	
	if(is.null(markers) && !is.null(covar.thresh)){

		for(i in 1:dim(covar.flags)[2]){
			is.covar <- which(t.stats[,i,] >= covar.thresh, arr.ind = TRUE)
			if(length(is.covar) > 0){
				covar.flags[is.covar[,1],i,is.covar[,2]] <- 1
				}
			}	
	
		data.obj$covar.flags <- covar.flags
		return(data.obj)
		}
	
	
	if(is.null(pheno)){
		pheno <- colnames(pheno.mat)
		}
	if(is.null(allele)){
		allele <- geno.names[[2]]
		}
	
	if(is.character(markers[1])){
		row.locale <- which(data.obj$geno.names[[3]] %in% markers)
		}else{
		row.locale <- markers
		}

	if(length(row.locale) < length(markers)){
		if(is.character(markers)[1]){
			didnt.find <- setdiff(markers, rownames(covar.flags)[row.locale])
			}else{
			didnt.find <- setdiff(markers, c(1:dim(covar.flags)[1])[row.locale])		
			}
		cat("\nI couldn't find the following markers:\n")
		cat(didnt.find, sep = "\n")
		return(data.obj)
		}


	if(is.character(pheno)[1]){
		pheno.locale <- which(colnames(covar.flags) %in% pheno)		
		}else{
		pheno.locale <- 1:dim(covar.flags)[2]
		}

	if(length(pheno.locale) < length(pheno)){
		if(is.character(pheno)[1]){
			didnt.find <- setdiff(pheno, colnames(covar.flags)[pheno.locale])
			}else{
			didnt.find <- setdiff(pheno, c(1:dim(covar.flags)[2])[pheno.locale])	
			}
		cat("\nI couldn't find the following phenotypes:\n")
		cat(didnt.find, sep = "\n")
		return(data.obj)	
		}

	if(is.character(allele)[1]){
		allele.locale <- which(dimnames(covar.flags)[[3]] %in% allele)
		}else{
		allele.locale <- 1:dim(covar.flags)[3]
		}

	if(length(allele.locale) < length(allele)){
		if(is.character(allele)[1]){
			didnt.find <- setdiff(allele, dimnames(covar.flags)[[3]][allele.locale])
			}else{
			didnt.find <- setdiff(allele, c(1:dim(covar.flags)[3])[allele.locale])	
			}
		cat("\nI couldn't find the following alleles:\n")
		cat(didnt.find, sep = "\n")
		return(data.obj)	
		}


	if(is.covar){
		covar.flags[row.locale, pheno.locale, allele.locale] <- 1
		}else{
		covar.flags[row.locale, pheno.locale, allele.locale] <- 0	
		}
		

	data.obj$covar.flags <- covar.flags
	
	pair.covar <- data.obj$covar.for.pairscan
	
	if(!is.null(pair.covar)){
		markers.in.pair <- which(rownames(covar.flags) %in% rownames(pair.covar))
		new.pairscan.covar <- covar.flags[markers.in.pair,,]
		data.obj$covar.for.pairscan <- new.pairscan.covar
		}

	return(data.obj)

}