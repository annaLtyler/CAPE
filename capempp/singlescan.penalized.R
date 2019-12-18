#This script performs a 1D scan of the data
#and can accept names to use as covariates.
#covariates must be encoded as markers (see create.covar.R)

#This script first calls the permutation script 
#genome.wide.significance.R to calculate the significance 
#threshold. It then uses this threshold to determine which 
#markers will be used as covariates for each phenotype

#The script creates a table of effect sizes for each marker
#for each phenotype, and a table of covariate flags indicating
#which markers will be used as covariates for each phenotoype
#These are both added to the data object

#The user has the choice to scan the eigentraits (default)
#or the original phenotypes. 

#ref.allele is used in a multi-allele scan. The ref.allele is
#removed from the regression, and all reported effects are in
#referece to the ref.allele.
#Threshold.fun and threshold.param work together to threshold
#the "significant" loci. threshold.fun might be "t.stat.thresh.sd",
#for example, which would threshold the t statistics by standard
#deviation. The threshold.param that would go with this, is how
#many standard deviations away from the mean you'd like to set 
#as the cutoff.

#data.obj <- cross.obj; n.perm = 10; scan.what = "raw.traits"; ref.allele = "B"; covar = c("Sex", "diet"); pairscan.thresh.fun = "p.val.thresh"; pairscan.thresh.param = 0.05; covar.thresh.fun = "p.val.thresh"; covar.thresh.param = 0.01

singlescan.penalized <- function(data.obj, geno.obj, n.perm = 100, scan.what = c("eigentraits", "raw.traits"), ref.allele = NULL, covar = NULL, alpha = c(0.01, 0.05), regression.type = c("standard", "ridge", "lasso", "elasticNet"), verbose = FALSE) {

	require(glmnet)
	
	if(length(regression.type) > 1){
		regression.type = "standard"
		}else{
		rtypes <- c("standard", "ridge", "lasso", "elasticNet")
		type.locale <- which(rtypes == regression.type)	
		if(length(type.locale) == 0){
			stop("regression.type must be one of 'standard', 'ridge', 'lasso', or 'elasticNet'")
			}	
		}

	check.underscore(data.obj)
	check.bad.markers(data.obj, geno.obj)

	geno.names <- data.obj$geno.names
		
	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentrais, otherwise, use raw phenotypes
	type.choice <- c(grep("eig", scan.what), grep("ET", scan.what), grep("et", scan.what))
	if(length(type.choice) > 0){
		pheno <- data.obj$ET
		}else{
			pheno <- data.obj$pheno
			}


	singlescan.obj <- vector(mode = "list", length = 7)
	names(singlescan.obj) <- c("alpha", "alpha.thresh", "ref.allele", "singlescan.effects", "singlescan.t.stats", "locus.p.vals", "locus.lod.scores")

	singlescan.obj$alpha <- alpha


	gene <- get.geno(data.obj, geno.obj)
	
	#Get the dimension names to minimize confusion	
	mouse.dim <- which(names(dimnames(gene)) == "mouse")
	locus.dim <- which(names(dimnames(gene)) == "locus")
	allele.dim <- which(names(dimnames(gene)) == "allele")

	if(regression.type != "standard"){
		if(verbose){cat("Flattening genotype matrix for", regression.type, "regression...\n")}
		flat.geno <- flatten.geno(gene, verbose = verbose)
		cat("\n")
		}
	
	n.phe <- dim(pheno)[2]
	n.gene <- dim(gene)[locus.dim]
	
	#first do the permutations to get the significance threshold
	#results will be compared to the significance threshold to 
	#determine which markers to use as covariates
	if(n.perm > 0 && regression.type == "standard"){
		if(verbose){cat("\n\nPerforming permutations to calculate significance threshold...\n")}
		singlescan.obj$alpha.thresh <- genome.wide.threshold.1D(data.obj, geno.obj, n.perm = n.perm, scan.what = scan.what, ref.allele = ref.allele, alpha = alpha, verbose = verbose)
		}else{
		if(regression.type == "standard"){
			message("Permutations are not being performed. No significance thresholds will be calculated.")
		}else{
			message("Permutations are not performed with ", regression.type, " regression.")
			}
		}

	#if there are covariates specified, pull these out.
	#covariates must be coded as markers and contained in the
	#genotype matrix
	if(!is.null(covar)){
		covar.vals <- data.obj$covar.table
		covar.table <- array(NA, dim = c(length(geno.names[[mouse.dim]]), length(geno.names[[allele.dim]]), length(covar)))
		for(i in 1:length(covar)){
			covar.table[,1:dim(covar.table)[2],i] <- covar.vals[,i]
			}
		}else{
			covar.table <- NULL
			}

	#check for a reference allele, pull it out of the 
	#allele names here and add it to the data object
	if(length(ref.allele) != 1){ #add a check for the reference allele
		stop("You must specify one reference allele")
		}
	ref.col <- which(dimnames(gene)[[allele.dim]] == ref.allele)
	if(length(ref.col) == 0){
		stop("I can't find reference allele: ", ref.allele)
		}
	new.allele.names <- dimnames(gene)[[allele.dim]][-ref.col]
	singlescan.obj$ref.allele <- ref.allele

	
	#=====================================================================
	#functions used to do regressions and gather stats
	#=====================================================================	


		#This function takes the results from get.stats.multiallele
		#and parses them into the final arrays
		add.results.to.array <- function(result.array, results.list, stat.name){
			row.num <- which(rownames(results.list[[1]][[1]]) == stat.name)
			#find the next spot to put data
			#Each successive phenotype is stored
			#in the 2nd dimension of the array
			next.spot.locale <- min(which(is.na(result.array[1,,1])))
			result.mat <- t(sapply(results.list, function(x) as.vector(x[[1]][row.num,])))
			result.array[,next.spot.locale,] <- t(sapply(results.list, function(x) as.vector(x[[1]][row.num,])))
			return(result.array)	
			}
			
		add.flat.results.to.array <- function(result.array, model, pheno.num){
			betas <- as.vector(model[[2]]$beta)
			num.markers <- dim(result.array)[1]
			num.alleles <- dim(result.array)[3]
			start.pos <- 1
			for(i in 1:num.markers){
				locus.results <- betas[start.pos:(start.pos+num.alleles-1)]
				result.array[i,pheno.num,] <- locus.results
				start.pos = start.pos + num.alleles
				}
			return(result.array)
			}
		#=====================================================================
		#end of internal functions
		#=====================================================================
	

		#=====================================================================
		#begin code for multi-allelic cross
		#=====================================================================
		#In the multi-allele case, we want to collect
		#three 3D arrays each of num.marker by num.pheno by num.allele:
		#array of t statistics (for plotting p vals of regressions)
		#array of effects (betas) (for effect plots)
		#array of covar flags (for use in pair.scan)

		t.stat.array <- array(dim = c(dim(gene)[[locus.dim]], dim(pheno)[2], (dim(gene)[[allele.dim]]-1)))
		effect.array <- array(dim = c(dim(gene)[[locus.dim]], dim(pheno)[2], (dim(gene)[[allele.dim]]-1)))
		dimnames(t.stat.array) <- dimnames(effect.array) <- list(dimnames(gene)[[locus.dim]], dimnames(pheno)[[2]], new.allele.names)
		
		#make arrays to hold locus-by-locus stats
		locus.lod.scores <- matrix(1, ncol = n.phe, nrow = dim(gene)[[locus.dim]])
		locus.p.values <- matrix(1, ncol = n.phe, nrow = dim(gene)[[locus.dim]])
		rownames(locus.lod.scores) <- rownames(locus.p.values) <- dimnames(gene)[[locus.dim]]
		colnames(locus.lod.scores) <- colnames(locus.p.values) <- colnames(pheno)

		for(i in 1:n.phe){
			if(verbose){cat("\nScanning trait:", colnames(pheno)[i], "\n")}
			#take out the response variable
			phenotype <- pheno[,i]
				
			#apply the modeling function to each marker column
			if(regression.type == "standard"){
			
				results.list <- apply(gene, locus.dim, function(x) get.stats.multiallele(phenotype, x, covar.table))
				#fill each result array with the returned values
				t.stat.array <- add.results.to.array(t.stat.array, results.list, "t.stat")
				effect.array <- add.results.to.array(effect.array, results.list, "slope")
				locus.lod.scores[,i] <- unlist(lapply(results.list, function(x) x$lod))
				locus.p.values[,i] <- unlist(lapply(results.list, function(x) x$pval))
			
				}else{
			
				if(regression.type == "ridge"){num.alpha = 1; alpha = 0}
				if(regression.type == "lasso"){num.alpha = 1; alpha = 1}
				if(regression.type == "elasticNet"){num.alpha = 10; alpha  = NULL}
				model <- el.net(flat.geno, pheno[,i], num.alpha = num.alpha, alpha = alpha)
				effect.array <- add.flat.results.to.array(effect.array, model, i)
				t.stat.array <- effect.array
				}
			}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
			if(verbose){cat("\n")}
		
			singlescan.obj$singlescan.effects <- effect.array
			singlescan.obj$singlescan.t.stats <- t.stat.array
			singlescan.obj$locus.p.vals <- locus.p.values
			singlescan.obj$locus.lod.scores <- locus.lod.scores
		

		return(singlescan.obj)
	
	}

