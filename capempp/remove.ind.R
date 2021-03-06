#This function removes individuals

remove.ind <- function(data.obj, ind.to.remove = NULL, names.to.remove = NULL){

		ind.idx <- NULL

		if(!is.null(names.to.remove)){		
			ind.idx <- match(names.to.remove, rownames(data.obj$pheno))
			}
		if(!is.null(ind.to.remove)){
			ind.idx = unique(c(ind.idx, ind.to.remove))
			}

		if(length(ind.idx) > 0){
			#remove individuals from the phenotype matrix in data.obj
			#and geno.names
			data.obj$pheno <- data.obj$pheno[-ind.idx,,drop=FALSE]
			data.obj$geno.names[[1]] <- data.obj$geno.names[[1]][-ind.idx]

			#if covariates have already been assigned, remove individuals
			#from these tables as well.
			if(!is.null(data.obj$p.covar.table)){
				data.obj$p.covar.table <- data.obj$p.covar.table[-ind.idx,,drop=FALSE]
				}
			if(!is.null(data.obj$g.covar.table)){
				data.obj$g.covar.table <- data.obj$g.covar.table[-ind.idx,,drop=FALSE]
				}
			if(!is.null(data.obj$raw.pheno)){
				data.obj$raw.pheno <- data.obj$raw.pheno[-ind.idx,,drop=FALSE]
				}
			if(!is.null(data.obj$ET)){
				warning("get.eigentraits needs to be re-run because individuals were removed.\n")	
				}
			}
		return(data.obj)
		}