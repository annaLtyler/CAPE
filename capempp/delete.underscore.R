#This is an internal function that checks for underscores
#in marker names and deletes them. We use it in a couple places
#to make 

delete.underscore <- function(data.obj, geno.obj = NULL){
	
	
	geno <- get.geno(data.obj, geno.obj)
	
	marker.names <- data.obj$geno.names[[3]]
	under.locale <- grep("_", marker.names)
	
	if(length(under.locale) > 0){
		bad.names <- marker.names[under.locale]
		new.names <- unlist(lapply(strsplit(bad.names, "_"), function(x) paste(x[1:length(x)], collapse = "")))
		
		data.obj$geno.names[[3]][under.locale] <- new.names
		dimnames(geno)[[3]][under.locale] <- new.names
		message("Removing underscores from marker names\n")
		# cat(bad.names, sep = "\n")
		}	
	
	results <- list(data.obj, geno)
	names(results) <- c("data.obj", "geno.obj")
	
	return(results)

	}

	