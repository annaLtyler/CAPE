#This function takes in a singlescan object
#and selects markers for the pairscan based on
#the effect sizes of the markers originally selected
#for the pairscan
#Currently, this function is only used in generating a null
#distribution, not in selecting markers initially. So it is
#essentially not used anymore. At some point, I'd like to use
#is for initial marker selection as well.
#!!!This function does not yet handle multi-parent crosses!!!

select.markers.for.pairscan.dist <- function(data.obj, singlescan.obj, geno.obj, verbose = TRUE){

	require(abind)
	
	if(is.null(data.obj$marker.selection.method)){
		data.obj$marker.selection.method <- "effects.dist"
		}
	#===============================================================
	# get the distribution of effect sizes across the traits
	# we want to match this distribution
	#===============================================================
	geno.for.pairscan <- data.obj$geno.for.pairscan
	split.markers <- strsplit(colnames(geno.for.pairscan), "_")
	just.markers <- unlist(lapply(split.markers, function(x) x[1]))
	just.alleles <- unlist(lapply(split.markers, function(x) x[2]))
	
	if(class(singlescan.obj) == "list"){ 
		results <- abs(singlescan.obj$singlescan.t.stats) #an actual singlescan object
		}else{
		results <- abs(singlescan.obj) #a singlescan matrix for calculating pairscan null distribution
		}

	filtered.results <- results
	
	covar.info <- get.covar(data.obj)
	results.no.covar <- results[which(!rownames(results) %in% covar.info$covar.names),,,drop=FALSE]

	selected.effects <- t(mapply(function(x, y) results.no.covar[x,,y,drop=FALSE], just.markers, just.alleles))
	if(class(singlescan.obj) == "list"){
		colnames(selected.effects) <- dimnames(singlescan.obj$singlescan.effects)[[2]]
		}else{
		colnames(selected.effects) <- dimnames(singlescan.obj)[[2]]	
		}

	orig.effects <- selected.effects	


	#===============================================================
	# find all effects in the current single scan
	#===============================================================
	
	if(class(singlescan.obj) == "list"){
		ref.allele <- singlescan.obj$ref.allele
		data.obj$ref.allele <- ref.allele
		}else{
		ref.allele <- data.obj$ref.allele
		}

	if(is.null(ref.allele)){
		allele.text <- paste(alleles, collapse = ", ")
		ref.allele <- readline(prompt = paste("Which allele do you want to use as the reference?\n", allele.text, "\n"))
		data.obj$ref.allele <- ref.allele
		}

	selected.markers <- vector(mode = "list", length = ncol(orig.effects))
	# selected.effects <- vector(mode = "list", length = ncol(orig.effects))
	for(i in 1:length(selected.markers)){
		selected.markers[[i]] <- unlist(lapply(orig.effects[,i], function(x) rownames(results.no.covar)[get.nearest.pt(results.no.covar[,i,], x)]))
		# selected.effects[[i]] <- unlist(lapply(orig.effects[,i], function(x) results.no.covar[get.nearest.pt(results.no.covar[,i,], x),i,]))		
		}
	u_selected <- unique(unlist(selected.markers))		
	
	# plot(density(orig.effects))
	# points(density(unlist(selected.effects)), col = "red", type = "l")
	
	geno <- get.geno(data.obj, geno.obj)
	marker.locale <- sort(match(u_selected, dimnames(geno)[[3]]))
	
	geno.for.pairscan <- geno[,ref.allele,marker.locale]
	
	#===============================================================
	# check the final list for linear independence
	#===============================================================

	if(verbose){cat("Checking for linear independence...\n")}
	data.obj$geno.for.pairscan <- geno.for.pairscan
	geno.ind <- get.linearly.independent(data.obj)
	
	rownames(geno.ind$independent.markers) <- rownames(data.obj$pheno)
	data.obj$geno.for.pairscan <- geno.ind$independent.markers
	
	if(verbose){
		cat(length(geno.ind[[2]]), "allele(s) rejected.\n")
		cat("Final alleles selected:", "\t", ncol(geno.ind$independent.markers), "\n")
		}
	data.obj$marker.selection.method = "effects.dist"				
	return(data.obj)
	
	
}