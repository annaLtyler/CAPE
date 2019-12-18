adjust.p.directed.influences <- function(data.obj, adjust.method = c("holm", "fdr")){

	dir.inf <- data.obj$max.var.to.pheno.influence
	full.length <- sum(unlist(lapply(dir.inf, nrow)))
	chunk.p <- chunkV(1:full.length, length(dir.inf))
	
	all.p <- unlist(lapply(dir.inf, function(x) as.numeric(x[,"emp.p"])))
	adj.p <- p.adjust(all.p, adjust.method)
	
	for(i in 1:length(dir.inf)){
		dir.inf[[i]][,"p.adjusted"] <- adj.p[chunk.p[[i]]]
		}
	data.obj$max.var.to.pheno.influence <- dir.inf
	return(data.obj)
	}
