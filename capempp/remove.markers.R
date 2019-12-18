#This function removes markers


remove.markers <- function(data.obj, markers.to.remove){

	marker.idx <- match(markers.to.remove, data.obj$geno.names[[3]])
	if(any(is.na(marker.idx))){marker.idx <- markers.to.remove}
		
	#take out the markers from the data.obj meta data	
	data.obj$chromosome <- data.obj$chromosome[-marker.idx]
	data.obj$marker.num <- data.obj$marker.num[-marker.idx]
	data.obj$marker.location <- data.obj$marker.location[-marker.idx]
	data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][-marker.idx]
	
	return(data.obj)
	}
	
