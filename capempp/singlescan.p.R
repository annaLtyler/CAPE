#This function gets actual p values for 
#a singlescan

singlescan.p <- function(data.obj, geno.obj, scan.what = c("ET", "raw"), n.cores = 4, run.parallel = TRUE){
	
	
	gene <- get.geno(data.obj, geno.obj)
	pheno <- get.pheno(data.obj, scan.what)			
	n.phe = dim(pheno)[2]
	chr.which <- unique(data.obj$chromosome)
	
	#get the covariates and assign the variables
	#to the local environment
	covar.info <- get.covar(data.obj)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}
		
	n.covar <- length(covar.names)

	flat.geno <- gene[,1,]


	geno.chunks <- chunkV(1:ncol(flat.geno), n.cores)
	
	get.p <- function(g.idx, pheno.v){
		if(is.null(covar.table)){
			all.models <- lapply(g.idx, function(x) lm(pheno.v~flat.geno[,x]))
			}else{
			all.models <- lapply(g.idx, function(x) lm(pheno.v~covar.table+flat.geno[,x]))	
			}
		p.col <- lapply(all.models, function(x) anova(x)$"Pr(>F)")
		p.val <- unlist(lapply(p.col, function(x) tail(x[which(!is.na(x))], 1)))
		return(p.val)
		}

	all.p.mat <- matrix(NA, ncol = ncol(pheno), nrow = ncol(flat.geno))

	for(p in 1:ncol(pheno)){
		cat('Scanning', colnames(pheno)[p], "\n")

	  if (run.parallel) {

	    cl <- makeCluster(n.cores)
	    registerDoParallel(cl)
	    all.p <- foreach(x = 1:length(geno.chunks), .export = c("covar.table", "flat.geno")) %dopar% {
	      get.p(geno.chunks[[x]], pheno[,p])
	    }
	    stopCluster(cl)

	  } else {

	    all.p <- c()
	    index <- 1:length(geno.chunks)
	    for (x in index) {
	      all.p <- cbind(all.p, get.p(geno.chunks[[x]], pheno[,p]))
	    }

	  }
	
		all.p.mat[,p] <- unlist(all.p)
		
		}

	return(all.p.mat)
	
	
}