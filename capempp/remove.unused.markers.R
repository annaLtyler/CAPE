#This function takes out markers on the sex chromosomes
#as well as invariant markers.

remove.unused.markers <- function(data.obj, geno.obj){
	
	#we do not scan markers on the sex chromosomes
	#take these out here.
	x.locale <- grep("X", data.obj$chromosome, ignore.case = TRUE)
	if(length(x.locale) > 0){
		cat("Removing markers on the X chromosome\n")
		data.obj$chromosome <- data.obj$chromosome[-x.locale]
		data.obj$marker.location <- data.obj$marker.location[-x.locale]
		data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][-x.locale]
		data.obj$marker.num <- data.obj$marker.num[-x.locale]		
		}
		
	y.locale <- grep("Y", data.obj$chromosome, ignore.case = TRUE)
	if(length(y.locale) > 0){
		cat("Removing markers on the Y chromosome\n")
		data.obj$chromosome <- data.obj$chromosome[-y.locale]
		data.obj$marker.location <- data.obj$marker.location[-y.locale]
		data.obj$marker.num <- data.obj$marker.num[-y.locale]	
		data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][-y.locale]
		}

	m.locale <- grep("M", data.obj$chromosome, ignore.case = TRUE)
	if(length(m.locale) > 0){
		cat("Removing markers on the Mitochondrial chromosome\n")
		data.obj$chromosome <- data.obj$chromosome[-m.locale]
		data.obj$marker.location <- data.obj$marker.location[-m.locale]
		data.obj$marker.num <- data.obj$marker.num[-m.locale]	
		data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][-m.locale]
		}

	
	#take out markers with missing data
	gene <- get.geno(data.obj, geno.obj)
	num.allele <- apply(gene, 3, function(x) sum(x, na.rm = TRUE))
	mono.allele <- which(num.allele == 0)
	if(length(mono.allele) > 0){
		message(paste("\nRemoving", length(mono.allele), "markers with no data:"))
		cat(data.obj$geno.names[[3]][mono.allele], sep = ", ")
		data.obj$chromosome <- data.obj$chromosome[-mono.allele]
		data.obj$marker.location <- data.obj$marker.location[-mono.allele]
		data.obj$marker.num <- data.obj$marker.num[-mono.allele]
		data.obj$geno.names[[3]] <- data.obj$geno.names[[3]][-mono.allele]
		}

	return(data.obj)

	}
