#This function finds a significant marker within
#a given block 


find.sig.marker <- function(data.obj, block.name, pheno.name){
	
	block.locale <- which(names(data.obj$linkage.blocks.collapsed) == block.name)
	block.markers <- data.obj$linkage.blocks.collapsed[[block.locale]]
	
	full.block.names <- names(data.obj$linkage.blocks.full[match(block.markers, unlist(data.obj$linkage.blocks.full))])
	
	full.net <- data.obj$full.net
	pheno <- full.net[,(nrow(full.net)+1):ncol(full.net)]
	
	pheno.locale <- which(colnames(pheno) == pheno.name)
	
	main.effects <- which(pheno[full.block.names,pheno.locale] != 0)
	
	block.markers <- unlist(data.obj$linkage.blocks.full[match(names(main.effects), names(data.obj$linkage.blocks.full))])
	
	return(block.markers)

	}