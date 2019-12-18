#This script performs a single 1D scan
#it is currently used in generating the
#2D null distribution
#covar.mat is a matrix of the covariates

one.singlescanDO <- function(phenotype.vector, genotype.mat, ref.allele = "B", model.family, covar.vector = NULL, run.parallel = TRUE, n.cores = 4){

	if(!run.parallel){n.cores = 1}

	gene <- genotype.mat
		
	#Get the dimension names to minimize confusion	
	mouse.dim <- which(names(dimnames(gene)) == "mouse")
	locus.dim <- which(names(dimnames(gene)) == "locus")
	allele.dim <- which(names(dimnames(gene)) == "allele")


	#if there are covariates specified, pull these out.
	#covariates must be coded as markers and contained in the
	#genotype matrix
	if(!is.null(covar.vector)){
		covar <- names(covar.vector[which(covar.vector == 1)])
		if(length(covar) > 0){
			covar.loc <- get.col.num(gene, covar, locus.dim)
			covar.table <- array(NA, dim = c(dim(gene)[[mouse.dim]], dim(gene)[[allele.dim]], length(covar)))
			covar.table[,,1:length(covar)] <- gene[,,covar.loc]
			}else{
				covar.table = NULL
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


		#=====================================================================
		#begin code for multi-allelic cross
		#=====================================================================
			
		#apply the modeling function to each marker column
		if (run.parallel) {

		  cl <- makeCluster(n.cores)
		  registerDoParallel(cl)
		  results <- foreach(m = 1:dim(gene)[[locus.dim]], .export = c("get.stats.multiallele", "check.geno")) %dopar% {
		    get.stats.multiallele(phenotype.vector, gene[,,m], covar.table, model.family, ref.col)
		  }
		  stopCluster(cl)

		} else {

		  results <- c()
		  index <- 1:dim(gene)[[locus.dim]]
		  for (m in index) {
		    results <- cbind(results, get.stats.multiallele(phenotype.vector, gene[,,m], covar.table, model.family, ref.col))
		  }

		}
		
		t.stat.array <- matrix(unlist(lapply(results, function(x) x[[1]]["t.stat",])), ncol = length(new.allele.names), byrow = TRUE)
		
		colnames(t.stat.array) <- new.allele.names
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
		
		return(t.stat.array)		
	}

