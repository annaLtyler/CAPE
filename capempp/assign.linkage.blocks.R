#This function assigns markers to linkage blocks
#based on peak assignments from diagonal means.

assign.linkage.blocks <- function(data.obj, chr, l.blocks){
	
	chr.locale <- which(data.obj$chromosome == chr)
	chr.markers <- data.obj$geno.names[[3]][chr.locale]
	
	start.block = 1
	block.locale <- sapply(1:max(l.blocks), function(x) which(l.blocks == x))
	
	linkage.blocks <- vector(mode = "list", length = max(l.blocks))

	#outline each block
	for(b in 1:length(block.locale)){		
		end.block <- start.block + (length(block.locale[[b]])/2)
		linkage.blocks[[b]] <- chr.markers[start.block:end.block]
		start.block <- end.block
		} 
	
	return(linkage.blocks)
	
	
}