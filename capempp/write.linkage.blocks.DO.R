#This function writes out the marker names in each linkage block

write.linkage.blocks.DO <- function(data.obj, filename = "Linkage.Blocks.txt", pos.filename = "Linkage.Blocks.Positions.txt", print.internal.markers = TRUE, sig.figs = 3){
	
	blocks <- data.obj$linkage.blocks.collapsed
	
	if(is.null(blocks)){
		stop("linkage.blocks() or get.network() must be run before writing the linkage blocks.\n")
		}
	
	init.marker <- sapply(blocks, function(x) x[1])
	end.marker <- sapply(blocks, function(x) tail(x, 1))
	
	init.pos <- get.marker.location(data.obj, sapply(strsplit(init.marker, "_"), function(x) x[1]))
	end.pos <- get.marker.location(data.obj, sapply(strsplit(end.marker, "_"), function(x) x[1]))
	
	max.len <- max(sapply(blocks, length))
	
	name.mat <- matrix(NA, ncol = max.len, nrow = length(blocks))
	rownames(name.mat) <- names(blocks)

	pos.mat <- cbind(init.pos, end.pos)
	rownames(pos.mat) <- names(blocks)

	marker.mat <- list2Matrix(blocks)	
	
	
	write.table(marker.mat, file = filename, sep = "\t", quote = FALSE, na = "", col.names = FALSE)
	write.table(pos.mat, file = pos.filename, sep = "\t", quote = FALSE, na = "", col.names = FALSE)
}