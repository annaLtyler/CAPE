#This function runs a cape analysis from a parameter file
#The input is the initial cape data object and a parameter
#file name
#This function assumes you already have all required libraries
#and functions loaded 
#if using gene.based marker selection, there must be a file in
#the working directory called gene.list.txt with an ordered list
#of genes in a column
#kinship.type can be either "overall" or "LTCO"
# parameter.file = "cape.parameters.txt"; p.or.q = 0.05; results.file = "cross.RData"; n.cores = 4; run.singlescan = TRUE; run.pairscan = TRUE; error.prop.coef = TRUE; error.prop.perm = TRUE; initialize.only = FALSE; verbose = TRUE; run.parallel = TRUE
#setwd("~/Documents/Data/Scleroderma/Results/Dominant_lung_pheno")

run.cape <- function(data.obj, geno.obj, kin.obj = NULL, parameter.file = "cape.parameters.txt", 
p.or.q = 0.05, path = ".", results.file = "cross.RData", n.cores = 4, run.singlescan = TRUE, 
run.pairscan = TRUE, error.prop.coef = TRUE, error.prop.perm = TRUE, initialize.only = FALSE, 
verbose = TRUE, run.parallel = TRUE, save.results = TRUE, plot.results = TRUE){
 
 
  results.file <- file.path(path, results.file)

	data.obj <- compare.markers(data.obj, geno.obj)
		
	results.base.name <- basename(gsub(".RData", "", results.file))

	parameter.file.name <- file.path(path, parameter.file)
	parameter.table <- readParameters(parameter.file.name)
	
	for(i in 1:nrow(parameter.table)){
		if(is.na(parameter.table[i,1])){
			assign(rownames(parameter.table)[i], NULL)	
			}else{
			vals <- strsplit(parameter.table[i,1], " ")[[1]]
			if(!is.na(suppressWarnings(as.numeric(vals[1])))){
				assign(rownames(parameter.table)[i], as.numeric(vals))
				}else{
				assign(rownames(parameter.table)[i], vals)
				}
			}
		}
	data.obj$ref.allele <- ref.allele
	
	#===============================================================
	# Determine whether the kinship object has already been calculated
	# If so, read it in. If not, calculate the kinship object and
	# write out the file
	#===============================================================

	if(as.logical(use.kinship)){
		kin.file <- file.path(path, paste0("kinship_", results.base.name, ".RData"))
		if(!file.exists(kin.file)){
			kin.obj <- Kinship(data.obj, geno.obj, type = kinship.type, pop = pop)
			if(save.results){saveRDS(kin.obj, kin.file)}
		}
	}
						
		#===============================================================
		# We need a complete genotype matrix to calculate the kinship
		# adjusted genotypes later on. 
		# Check for missing values in the genotype matrix.
		# If there are missing values, impute them. 
		# Write out the imputed matrix, or read this in if it already
		# exists.
		#===============================================================
				
		#first check for imputed data. If it exists, read it in.
		imp.data.file <- file.path(path, paste0(results.base.name, "_data_imputed.RData"))
		imp.geno.file <- file.path(path, paste0(results.base.name, "_geno_imputed.RData"))
	
		#If it exists, read it in.
		if(as.logical(use.kinship)){
		if(file.exists(imp.geno.file)){
			data.obj <- readRDS(imp.data.file)
			geno.obj <- readRDS(imp.geno.file)		
		}else{ #otherwise, check the genotype object for missing data
			geno <- get.geno(data.obj, geno.obj)
			missing.vals <- which(is.na(geno))
			if(length(missing.vals) > 0){ #if there are missing values, impute them
				cat("There are missing values in geno.obj. Running impute.missing.geno...\n")
				geno.imp <- impute.missing.geno(data.obj, geno.obj)
				data.obj <- geno.imp$data.obj
				if(save.results){saveRDS(data.obj, imp.data.file)}
				geno.obj <- geno.imp$geno.obj
				if(save.results){saveRDS(geno.obj, imp.geno.file)}
				} # end case for when there are missing genotype values
			} #end case for when imputed files do not exist		
		#recalculate kinship with imputed data
		kin.obj <- Kinship(data.obj, geno.obj, type = kinship.type, pop = pop)
		if(save.results){saveRDS(kin.obj, kin.file)}
		}
		
		
		
		
	if(save.results){saveRDS(data.obj, results.file)}

	if(any(!run.singlescan, !run.pairscan, !error.prop.coef, !error.prop.perm)){
		data.obj <- readRDS(results.file)
		}else{

		if(verbose){cat("Removing unused markers...\n")}
		data.obj <- remove.unused.markers(data.obj, geno.obj)
		combined.data.obj <- delete.underscore(data.obj, geno.obj)
		
		data.obj <- combined.data.obj$data.obj
		geno.obj <- combined.data.obj$geno.obj
					
		if(!is.null(covariates)){
			data.obj <- pheno2covar(data.obj, covariates)
			}
		if(!is.null(marker.covariates)){
			data.obj <- marker2covar(data.obj, geno.obj, markers = marker.covariates)
			}
		
		data.obj <- select.pheno(data.obj, pheno.which = traits)	
	
		if(length(grep("e", scan.what, ignore.case = TRUE)) > 0){
			data.obj <- get.eigentraits(data.obj, scale.pheno = as.logical(traits.scaled), 
			normalize.pheno = as.logical(traits.normalized))
	
			if(plot.results){
				svd.file <- file.path(path, "svd.jpg")
				jpeg(svd.file, res = 300, width = 7, height = 7, units = "in")
				plotSVD(data.obj, orientation = "vertical", show.var.accounted = TRUE)
				dev.off()
				}
			
			data.obj <- select.eigentraits(data.obj, traits.which = eig.which)

			if(use.kinship){
				#if individuals were deleted from the phenotype matrix, delete these
				#from the kinship object too
				kin.obj <- remove.kin.ind(data.obj, kin.obj)
				}
			}
		
		if(save.results){saveRDS(data.obj, results.file)}
		}
	
	if(initialize.only){
		return(data.obj)
		}

	#===============================================================
	# run singlescan
	#===============================================================
	if(run.singlescan){
	singlescan.results.file <- file.path(path, paste0(results.base.name, ".singlescan.RData"))

		singlescan.obj <- singlescan(data.obj, geno.obj, kin.obj = kin.obj, n.perm = singlescan.perm, 
		ref.allele = ref.allele, alpha = c(0.01, 0.05), scan.what = scan.what, verbose = verbose, 
		run.parallel = run.parallel, n.cores = n.cores, model.family = "gaussian", overwrite.alert = FALSE)
		if(save.results){saveRDS(singlescan.obj, singlescan.results.file)}
		

		if(plot.results){
			for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
				fig.file <- file.path(path, paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], 
				".Standardized.jpg"))
				jpeg(fig.file, width = 20, height = 6, units = "in", res = 300)
				plotSinglescan(data.obj, singlescan.obj = singlescan.obj, standardized = TRUE, 
				allele.labels = NULL, alpha = c(0.05, 0.01), include.covars = TRUE, line.type = "l", 
				pch = 16, cex = 0.5, lwd = 3, traits = colnames(singlescan.obj$singlescan.effects)[ph])
				dev.off()
				}
	
			for(ph in 1:ncol(singlescan.obj$singlescan.effects)){
				fig.file <- file.path(path, paste0("Singlescan.", colnames(singlescan.obj$singlescan.effects)[ph], 
				".Effects.jpg"))
				jpeg(fig.file, width = 20, height = 6, units = "in", res = 300)			
				plotSinglescan(data.obj, singlescan.obj = singlescan.obj, standardized = FALSE, allele.labels = NULL, 
				alpha = c(0.05, 0.01), include.covars = TRUE, line.type = "l", pch = 16, cex = 0.5, lwd = 3, 
				traits = colnames(singlescan.obj$singlescan.effects)[ph])
				dev.off()
				}
			}


		}else{
		singlescan.obj <- readRDS(singlescan.results.file)
		}
			
	
	#===============================================================
	# run pairscan
	#===============================================================
	
	pairscan.file <- file.path(path, paste0(results.base.name, ".pairscan.RData"))
	if(run.pairscan){
		if(marker.selection.method == "top.effects"){
			data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj, 
			num.alleles = num.alleles.in.pairscan, peak.density = peak.density, verbose = verbose, 
			plot.peaks = FALSE)
			}
			
		if(marker.selection.method == "from.list"){
			snp.file <- file.path(path, SNPfile)
			specific.markers <- read.table(snp.file, sep = "\t", stringsAsFactors = FALSE)
			data.obj <- select.markers.for.pairscan(data.obj, singlescan.obj, geno.obj, 
			specific.markers = specific.markers[,1], verbose = verbose, plot.peaks = FALSE)
			}
		
			
		if(marker.selection.method == "uniform"){
			data.obj <- select.markers.for.pairscan.uniform(data.obj, geno.obj, ref.allele = ref.allele, 
			required.markers = NULL, num.alleles = num.alleles.in.pairscan, verbose = verbose)
			}
	
	
		if(marker.selection.method == "by.gene"){
			gene.file <- file.path(path, "gene.list.txt")
			gene.list.mat <- read.table(gene.file, sep = "\t", stringsAsFactors = FALSE)		
			gene.list <- gene.list.mat[,1]
			data.obj <- select.markers.for.pairscan.by.gene(data.obj, geno.obj, ref.allele = ref.allele, 
			gene.list = gene.list, num.snps = num.alleles.in.pairscan, bp.buffer = bp.buffer, 
			organism = organism)
        } else {
			gene.list <- NULL
        }
	
		if(save.results){saveRDS(data.obj, results.file)}
	
		pairscan.obj <- pairscan(data.obj, geno.obj, scan.what = scan.what, 
		pairscan.null.size = pairscan.null.size, min.per.genotype = min.per.geno, 
		max.pair.cor = max.pair.cor, verbose = verbose, num.pairs.limit = Inf, 
		overwrite.alert = FALSE, run.parallel = run.parallel, n.cores = n.cores, 
		gene.list = gene.list, kin.obj = kin.obj)
		
		if(save.results){saveRDS(pairscan.obj, pairscan.file)}
					
		if(save.results){saveRDS(data.obj, results.file)}
    } else {
		pairscan.obj <- readRDS(pairscan.file)
    }
	
		
	#===============================================================
	# run reprametrization
	#===============================================================

	if(error.prop.coef){
		data.obj <- error.prop(data.obj, pairscan.obj, perm = FALSE, verbose = verbose, n.cores = n.cores, 
		run.parallel = run.parallel)
		if(save.results){saveRDS(data.obj, results.file)}
		}
	
	if(error.prop.perm){	
		data.obj <- error.prop(data.obj, pairscan.obj, perm = TRUE, verbose = verbose, n.cores = n.cores, 
		run.parallel = run.parallel)
		if(save.results){saveRDS(data.obj, results.file)}
		}
	
	data.obj <- calc.p(data.obj, pval.correction = pval.correction)
	
	if(length(grep("e", scan.what, ignore.case = TRUE)) > 0){
		transform.to.phenospace <- TRUE
		}else{
		transform.to.phenospace <- FALSE	
		}
		
		
	if(save.results){save.permutations = TRUE}else{save.permutations = FALSE}
	data.obj <- direct.influence(data.obj, pairscan.obj, transform.to.phenospace = transform.to.phenospace, 
	verbose = TRUE, pval.correction = pval.correction, save.permutations = save.permutations, 
	n.cores = n.cores, path = path)
	
	if(save.results){saveRDS(data.obj, results.file)}

	if(save.results){
		var.inf.file <- file.path(path, "Variant.Influences.csv")
		writeVariantInfluences(data.obj, p.or.q = p.or.q, filename = var.inf.file)
		}
	
	if(plot.results){
		fig.file <- file.path(path, "variant.influences.jpg")
		jpeg(fig.file, width = 10, height = 7, units = "in", res = 300)
		plotVariantInfluences(data.obj, p.or.q = p.or.q, standardize = FALSE, not.tested.col = "lightgray", 
		covar.width = 30, pheno.width = 30)
		dev.off()
		}
	
	data.obj <- get.network(data.obj, p.or.q = p.or.q, collapse.linked.markers = FALSE)
	data.obj <- get.network(data.obj, p.or.q = p.or.q, threshold.power = 1, collapse.linked.markers = TRUE, 
	plot.linkage.blocks = FALSE)
	if(save.results){saveRDS(data.obj, results.file)}
	
	if(plot.results){
		fig.file <- file.path(path, "Network.Circular.jpg")
		jpeg(fig.file, height = 7, width = 7, units = "in", res = 300)
		plotNetworkDO(data.obj, label.gap = 10, label.cex = 1.5, show.alleles = FALSE)
		dev.off()
		
		if(dim(geno.obj)[2] == 8){
			fig.file <- file.path(path, "Network.Circular.DO.jpg")
			jpeg(fig.file, height = 7, width = 7, units = "in", res = 300)
			plotNetworkDO(data.obj, label.gap = 10, label.cex = 1.5, show.alleles = TRUE)
			dev.off()		
			}	
		
		fig.file <- file.path(path, "Network.View.jpg")
		jpeg(fig.file, height = 7, width = 7, units = "in", res = 300)
		net.layout <- plotNetwork2(data.obj, zoom = 1.2, node.radius = 0.3, label.nodes = TRUE, 
		label.offset = 0.4, label.cex = 0.5, bg.col = "lightgray", arrow.length = 0.1, 
		layout.matrix = "layout_with_kk", legend.position = "topright", edge.lwd = 1, legend.radius = 2, 
		legend.cex = 0.7, xshift = -1)
		dev.off()
		}
	
	if(save.results){saveRDS(data.obj, results.file)}

	invisible(data.obj)

	
}
