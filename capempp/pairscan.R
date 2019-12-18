#This script performs the pairwise scan on all markers
#It takes in the data as a cross object.
#The user has the choice to scan the eigentraits (default)
#or the original phenotypes.
#This script also calls the function to do permutations
#on the 2D scan. It adds the genome-wide threshold for
#the 2D scan to the data object
#n.top.markers is used in generating the null. A permutation
#of the singlescan is run, and the n top markers are used
#in a permutation of the pairscan. if n.top.markers is NULL,
#it defaults to the number of markers in geno.for.pairscan
#marker.selection.method = c("top.effects", "uniform", "effects.dist", "by.gene" "by.list")

#ref.allele = "B"; pairscan.null.size = 1000; scan.what = "eigentraits"; max.pair.cor = 0.5; use.pairs.threshold = TRUE; verbose = TRUE; num.pairs.limit = 1e4; num.perm.limit = 1e7; use.kinship = FALSE; overwrite.alert = FALSE; kin.obj = NULL;min.per.genotype = NULL;run.parallel = TRUE


pairscan <- function(data.obj, geno.obj = NULL, ref.allele = "A", 
scan.what = c("eigentraits", "raw.traits"), pairscan.null.size = NULL, 
max.pair.cor = NULL, min.per.genotype = NULL, kin.obj = NULL, verbose = 
FALSE, num.pairs.limit = 1e6, num.perm.limit = 1e7, overwrite.alert = TRUE, 
run.parallel = TRUE, n.cores = 4, gene.list = NULL) {

	marker.selection.method = data.obj$marker.selection.method
	
	n.top.markers <- ncol(data.obj$geno.for.pairscan)

	if(!run.parallel){n.cores = 1}

	if(is.null(kin.obj)){use.kinship = FALSE}
	if(!is.null(kin.obj)){use.kinship = TRUE}
	data.obj$use.kinship <- use.kinship
	
	
	if(overwrite.alert){
		choice <- readline(prompt = "Please make sure you assign the output 
		of this function to a pairscan.obj, and NOT the data.obj. It will 
		overwrite the data.obj.\nDo you want to continue (y/n) ")
		if(choice == "n"){stop()}
		}


	pairscan.obj <- list()
	pairscan.obj$ref.allele <- ref.allele
	pairscan.obj$max.pair.cor <- max.pair.cor
	pairscan.obj$min.per.genotype <- min.per.genotype

	if(is.null(pairscan.null.size)){
		stop("The final size of the null distribution must be specified.")
		}
	
	pheno <- get.pheno(data.obj, scan.what)	
	num.pheno <- dim(pheno)[2]
	pheno.names <- colnames(pheno)

	covar.info <- get.covar(data.obj)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}
	#find the phenotypic covariates. These will
	#be tested separately, and not as part of a
	#chromosome
	num.covar <- length(covar.names)
	p.covar.locale <- which(covar.type == "p")
	num.p.covar <- length(p.covar.locale)
	

	if(is.null(data.obj$geno.for.pairscan)){
		stop("select.markers.for.pairscan() must be run before pairscan()")
		}
		
	
	#add the covariates (geno and pheno) 
	#to the genotype matrix so that we 
	#test all pairs
	gene <- get.geno.with.covar(data.obj, geno.obj, g.covar = TRUE, p.covar = TRUE, 
	for.pairscan = TRUE)	
	num.markers <- dim(gene)[2]
			
	#fill in a matrix to index the marker pairs
	if(verbose){cat("Getting marker pairs for pairscan...\n")}
	pared.marker.mat <- get.pairs.for.pairscan(gene, covar.names, max.pair.cor, 
	min.per.genotype, run.parallel = run.parallel, n.cores = n.cores, verbose = 
	verbose)

	num.pairs <- dim(pared.marker.mat)[1]
	
	if(num.pairs == 0){
		stop("There are no pairs to test. Try raising max.pair.cor or reducing 
		min.per.genotype.")
		}

	if(!is.null(num.pairs.limit) && num.pairs > num.pairs.limit){
		cat("\nThe number of pairs (",num.pairs,") exceeds ", num.pairs.limit, ".\n", sep = "")
		go.on <- readline(prompt = "Do you want to continue (y/n)?\n")
		if(length(grep("n", go.on))){
			message("Stopping pairwise scan...\n")
			return(pairscan.obj)
		}else{
			cat("Continuing pairwise scan...\n")
		}
	}

	cat("Performing pairwise tests...\n")
	#run one.pairscan for each phenotype with results in scanone.result
	if(!use.kinship){
		pairscan.results <- pairscan.noKin(data.obj, pheno.mat = pheno, 
		geno.mat = gene, covar.table = covar.table, paired.markers = pared.marker.mat, 
		n.perm = pairscan.null.size, verbose = verbose, n.cores = n.cores, 
		run.parallel = run.parallel)
		}else{
		pairscan.results <- pairscan.kin(data.obj, geno.obj = geno.obj, 
		scan.what = scan.what, marker.pairs = pared.marker.mat, kin.obj = kin.obj, 
		verbose = verbose, run.parallel = run.parallel, n.cores = n.cores)
		}	
		
		# print(str(pairscan.results))
		
	pairscan.obj$pairscan.results <- pairscan.results	

	if(pairscan.null.size > 0){	
		if(use.kinship){
		pairscan.perm <- pairscan.null.kin(data.obj, geno.obj, kin.obj, 
		scan.what = scan.what, ref.allele = ref.allele, pairscan.null.size = 
		pairscan.null.size, max.pair.cor = max.pair.cor, min.per.genotype, 
		verbose = verbose, marker.selection.method = marker.selection.method, 
		gene.list = gene.list, bp.buffer = bp.buffer, organism = organism)			
		}else{
		pairscan.perm <- pairscan.null(data.obj, geno.obj, scan.what = scan.what, 
		ref.allele = ref.allele, pairscan.null.size = pairscan.null.size, 
		max.pair.cor = max.pair.cor, min.per.genotype, verbose = verbose, 
		marker.selection.method = marker.selection.method, gene.list = gene.list, 
		bp.buffer = bp.buffer, organism = organism)
		}
	#add the results to the data object
	pairscan.obj$pairscan.perm <- pairscan.perm$pairscan.perm 
	#add the results to the data object
	pairscan.obj$pairs.tested.perm <- pairscan.perm$pairs.tested.perm 

		}
		
		
	return(pairscan.obj) #and return it
	
}
