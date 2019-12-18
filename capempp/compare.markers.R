compare.markers <- function(data.obj, geno.obj){	
	geno <- get.geno(data.obj, geno.obj)
	missing.markers <- setdiff(data.obj$geno.names[[3]], dimnames(geno)[[3]])
	if(length(missing.markers) > 0){
		message("Removing markers from data.obj that are not present in the geno.obj")
		data.obj <- remove.markers(data.obj, missing.markers)
		}
	return(data.obj)	
}
