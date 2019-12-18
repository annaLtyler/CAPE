#This function writes saves a bare-bones
#data object containing only the data used
#by the shiny app cape_interactions


write.data.for.shiny <- function(data.obj, name = "myData", path = "."){
		require(igraph)	
		
		shiny.obj <- list()
		shiny.obj$pheno <- data.obj$pheno
		shiny.obj$raw.pheno <- data.obj$raw.pheno
		shiny.obj$geno <- data.obj$geno.for.pairscan
	
		all.markers <- unlist(lapply(strsplit(colnames(data.obj$geno.for.pairscan), "_"), function(x) x[1]))
		marker.locale <- match(all.markers, data.obj$geno.names[[3]])
		
		shiny.obj$chromosome <- data.obj$chromosome[marker.locale]
		shiny.obj$marker.location <- data.obj$marker.location[marker.locale]
		shiny.obj$marker.num <- data.obj$marker.num[marker.locale]
		shiny.obj$geno.names <- data.obj$geno.names
		shiny.obj$geno.names[[3]] <- data.obj$geno.names[[3]][marker.locale]
		shiny.obj$var.to.var.p.val <- data.obj$var.to.var.p.val
		shiny.obj$max.var.to.pheno.influence <- data.obj$max.var.to.pheno.influence
		shiny.obj$p.covar <- data.obj$p.covar
		shiny.obj$p.covar.table <- data.obj$p.covar.table
		shiny.obj$g.covar <- data.obj$g.covar
		shiny.obj$g.covar.table <- data.obj$g.covar.table
	

	#========================================================================
	# code copied from plotVariantInfluences to make influence mat 
	# that can be thresholded later, since it takes a long time to 
	# fill it
	#========================================================================

	p.or.q = 1; standardize = TRUE
	geno.names <- data.obj$geno.names
	marker.names <- geno.names[[3]]

	var.influences <- data.obj$var.to.var.p.val
	pheno.inf <- data.obj$max.var.to.pheno.influence
	pheno.names <- names(data.obj$max.var.to.pheno.influence)
	num.pheno <- length(pheno.names)
	
	unique.markers <- unique(c(as.vector(var.influences[,"Source"]), as.vector(var.influences[,"Target"]), rownames(pheno.inf[[1]])))

		
	just.markers <- sapply(strsplit(unique.markers, "_"), function(x) x[[1]][1])
	unique.marker.locale <- match(just.markers, marker.names)
	sorted.markers <- unique.markers[order(unique.marker.locale)]

	#update the markers based on the covariate width
	orig.chromosomes <- data.obj$chromosome[sort(unique.marker.locale)]
	
	covar.info <- get.covar(data.obj)
	
					
	var.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	var.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	colnames(var.influence.mat) <- rownames(var.influence.mat) <- colnames(var.pval.mat) <- rownames(var.pval.mat) <- sorted.markers

	pheno.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	pheno.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	colnames(pheno.influence.mat) <- colnames(pheno.pval.mat) <- pheno.names
	rownames(pheno.influence.mat) <- rownames(pheno.pval.mat) <- sorted.markers
	
	
	#fill the variant-to-variant matrix with test statistics with sources in rows and targets in columns
	var.sig.col <- which(colnames(var.influences) == "p.adjusted")
	var.sig.net <- var.pval.net <- graph_from_edgelist(as.matrix(var.influences[,1:2]))
	 
	if(standardize){
		E(var.sig.net)$weight <- as.numeric(var.influences[,"Effect"])/as.numeric(var.influences[,"SE"])
		}else{
		E(var.sig.net)$weight <- as.numeric(var.influences[,"Effect"])
		}
	E(var.pval.net)$weight <- as.numeric(var.influences[,var.sig.col])
	
	var.pval.mat <- as.matrix(as_adjacency_matrix(var.pval.net, attr = "weight"))
	var.influence.mat <- as.matrix(as_adjacency_matrix(var.sig.net, attr = "weight"))	

	
	#fill the variant-to-phenotype matrix with test statistics 
	#(still with sources in rows and targets in columns)
	#use phenotypes or eigentraits based on user input
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.inf[[1]]))))
	for(i in 1:length(unique.markers)){
			for(j in 1:length(pheno.names)){
				marker.locale <- which(rownames(pheno.inf[[j]]) == unique.markers[i])
				if(length(marker.locale) > 0){
					if(standardize){	
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "t.stat"]
						}else{
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "coef"]
						}
				pheno.pval.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, pheno.sig.col]
				}else{
					pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- NA
					pheno.pval.mat[unique.markers[i], pheno.names[j]] <- NA
					}
			}
		}
		

	
	full.inf.mat <- cbind(var.influence.mat, pheno.influence.mat)
	full.pval.mat <- cbind(var.pval.mat, pheno.pval.mat) 

	full.inf.mat <- apply(full.inf.mat, 2, as.numeric)
	full.pval.mat <- apply(full.pval.mat, 2, as.numeric)
	
	rownames(full.inf.mat) <- sorted.markers
	colnames(full.inf.mat) <- c(sorted.markers, pheno.names)
	
	shiny.obj$var.inf.mat <- full.inf.mat
	shiny.obj$pval.mat <- full.pval.mat	


	#make a variant influences matrix that can be thresholded later
	
	saveRDS(shiny.obj, file = paste(path, "/", name, ".RData", sep = ""))
	
	invisible(shiny.obj)
	}