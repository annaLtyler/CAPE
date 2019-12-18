#This function writes the parameters 
#used in a cape analysis

writeParameters <- function(data.obj, singlescan.obj, pairscan.obj, filename = "cape.parameters.txt"){

	#================================================
	# internal functions
	#================================================
	write.header <- function(header.name){
		cat("#================================================\n")
		cat("#", header.name, "\n")
		cat("#================================================\n")
		}

	write.param <- function(param.list){
		for(i in 1:length(param.list)){
			if(!is.null(get(param.list[i]))){
				cat(param.list[i], "\t", get(param.list[i]), "\n")
				}
			}
		cat("\n")
		}
	
	#================================================
	# general parameters
	#================================================
	gen.param <- c("traits", "covariates", "marker.covariates", "traits.scaled", "traits.normalized", "scan.what", "eig.which", "use.kinship", "pval.correction")
		
	traits <- colnames(data.obj$pheno)
	covar <- data.obj$p.covar
	marker.covar <- data.obj$g.covar
	traits.scanned <- names(pairscan.obj$pairscan.results)
	eig.which <- 1:ncol(data.obj$ET)
	if(length(eig.which) > 0 || length(which(traits.scanned %in% phenotypes) > 0)){
		scan.what <- "eigentraits"
		}else{
		scan.what <- "raw.traits"
		}
	traits.scaled <- data.obj$traits.scaled
	traits.normalized <- data.obj$traits.normalized
	use.kinship <- data.obj$use.kinship
	pval.correction <- data.obj$pval.correction

	#================================================
	#single scan parameters
	#================================================
	single.param <- c("ref.allele", "singlescan.perm")
	
	ref.allele <- pairscan.obj$ref.allele
	singlescan.perm <- singlescan.obj$n.perm
	
	#================================================
	# marker selection
	#================================================
	marker.param <- c("marker.selection.method", "peak.density", "tolerance", "window.size", "num.alleles.in.pairscan", "bp.buffer", "organism")
	
	num.alleles.in.pairscan = ncol(data.obj$geno.for.pairscan)
	bp.buffer <- data.obj$bp.buffer
	marker.selection.method <- data.obj$marker.selection.method
	if(marker.selection.method == "top.effects"){
		peak.density = data.obj$peak.density
		tolerance = data.obj$tolerance
		window.size = data.obj$window.size
		}else{
		peak.density = NULL
		tolerance = NULL
		window.size = NULL
		}
	organism <- data.obj$organism

	#================================================
	# pair scan
	#================================================
	pair.param <- c("max.pair.cor", "min.per.geno", "pairscan.null.size")
	max.pair.cor = pairscan.obj$max.pair.cor
	min.per.geno <- pairscan.obj$min.per.geno
	pairscan.null.size = nrow(pairscan.obj$pairscan.perm[[1]][[1]])
	
	#================================================
	# write out all the parameter sets
	#================================================
	sink(filename)
	write.header("General Parameters")
	write.param(gen.param)

	write.header("Single Scan Parameters")
	write.param(single.param)

	write.header("Marker Selection Parameters")
	write.param(marker.param)

	write.header("Pairscan Parameters")
	write.param(pair.param)

	sink()


}