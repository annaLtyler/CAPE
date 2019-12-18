#This function takes in a genotype matrix and returns a list of all
#pairs that pass the pass.checks tests


get.pairs.for.pairscan <- function(gene, covar.names = NULL, max.pair.cor = NULL, min.per.genotype = NULL, run.parallel = TRUE, n.cores = 4, verbose = TRUE){

	if(!run.parallel){n.cores = 1}

	p = NULL #for appeasing R CMD check
	
	if(is.null(max.pair.cor) && is.null(min.per.genotype)){
		stop("One of max.pair.cor or min.per.genotype should be set.")
		}
	if(!is.null(max.pair.cor) && !is.null(min.per.genotype)){
		stop("Only one of max.pair.cor or min.per.genotype should be set.")
		}
	
	if(verbose){
		cat("\nChecking marker pairs for genotype representation...\n")
		}
	
	#covariates should be paired with every marker
	#so only check pairs for genetic markers, and make
	#a separate list of pairs for covariates with
	#markers and each other
	all.covar.pairs <- NULL
	if(!is.null(covar.names)){
		covar.locale <- match(covar.names, colnames(gene))
		covar.mat <- gene[,covar.locale,drop=FALSE]
		gene <- gene[,-covar.locale,drop=FALSE]
		covar.pairs <- pair.matrix(covar.names)
		covar.geno.pairs <- cbind(rep(colnames(gene), length(covar.names)), rep(covar.names, each = ncol(gene)))
		all.covar.pairs <- rbind(covar.geno.pairs, covar.pairs)
		}
	
	if(!is.null(max.pair.cor)){
		thresh.param <- max.pair.cor
		check.linkage <- function(m1,m2,thresh.param){
			pair.cor <- try(cor(m1, m2, use = "complete"), silent = TRUE)
			if(class(pair.cor) == "try-error" || pair.cor > max.pair.cor || is.na(pair.cor)) {
				return(FALSE) #pair failed check
				}else{
				return(TRUE) #pair passed check
				}
			}
	}
	
	
	if(!is.null(min.per.genotype)){
		thresh.param <- min.per.genotype
		check.linkage <- function(m1,m2,thresh.param){
			geno.table <- cbind.data.frame(as.factor(m1),as.factor(m2))
			colnames(geno.table) <- c("m1","m2")
			reps <- table(geno.table$m1,geno.table$m2)
			too.few <- which(reps < thresh.param)
			if(length(too.few) >= 1) {
				return(FALSE) #pair failed check
				}else{
				return(TRUE) #pair passed check
				}
			}		
		}

	
	all.pairs <- pair.matrix(1:dim(gene)[2])
	all.pair.names <- pair.matrix(colnames(gene))
	
	if(verbose){cat("There are", nrow(all.pairs), "possible marker pairs to test.\n")}
	
	check.one.pair <- function(pair.num){
		pair <- all.pairs[pair.num,]
		pass.checks <- check.linkage(m1 = gene[,pair[1]], m2 = gene[,pair[2]], thresh.param = thresh.param)
		return(pass.checks)
		}	
	
	check.multi.pairs <- function(pair.V){
		pair.checks <- unlist(lapply(pair.V, check.one.pair))
		# result <- cbind(pair.V, pair.checks)
		# return(result)
		return(pair.checks)
		}
	
	#chunk up the jobs based on how many cores we want to use
	pair.list <- chunkV(1:nrow(all.pairs), n.cores)
	#str(pair.list)
	
	if (run.parallel) {

	  cl <- makeCluster(n.cores)
	  registerDoParallel(cl)
	  good.pair.list <- foreach(p = 1:length(pair.list)) %dopar% {
	    check.multi.pairs(pair.list[[p]])
	  }
	  stopCluster(cl)

	} else {

	  good.pair.list <- c()
	  index <- 1:length(pair.list)
	  for (p in index) {
	    good.pair.list <- rbind(good.pair.list, check.multi.pairs(pair.list[[p]]))
	  }

	}
	
	testV <- unlist(good.pair.list)
	idxV <- unlist(pair.list)
	
	pairs.mat <- all.pair.names[idxV[which(testV)],,drop = FALSE]
	if(verbose){
		cat(dim(pairs.mat)[1], "marker pairs will be tested.\n")
		}
		
	colnames(pairs.mat) <- c("marker1", "marker2")
	rownames(pairs.mat) <- NULL
	pairs.mat <- rbind(pairs.mat, all.covar.pairs)
	
	if(verbose){
		cat(dim(pairs.mat)[1], "pairs including covariates will be tested.\n")
		}
	
	return(pairs.mat)
	
}