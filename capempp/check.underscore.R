#This is an internal function that checks for underscores
#in marker names and deletes them. We use it in a couple places
#to make 

check.underscore <- function(data.obj){
	
	marker.names <- data.obj$geno.names[[3]]
	under.locale <- grep("_", marker.names)
	if(length(under.locale) > 0){
		stop("Underscores have been detected in some marker names.\nPlease use delete.underscore() to remove these before proceeding.")
		}
	}

	