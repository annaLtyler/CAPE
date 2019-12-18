#This function writes out a cape cross
#in a csv format so you can modify a
#cross and write it out in a nice format

writePopulation <- function(data.obj, geno.obj, ref.allele = "A", na = NA, filename = "capeData.csv"){
	
  if(class(geno.obj)=="array"){geno.obj<-geno.obj}else{geno.obj<-geno.obj$geno}
  
	if(dim(geno.obj)[2] > 2){
		stop("writePopulation only works for two-parent populations.")
		}

	ref.allele <- ref.allele[1]
	alt.allele <- setdiff(colnames(geno.obj), ref.allele)
	geno <- get.geno(data.obj, geno.obj)
	geno.mat <- geno[,which(colnames(geno) == alt.allele),]
	pheno <- data.obj$pheno
	chr <- data.obj$chromosome
	loc <- data.obj$marker.location
	marker.names <- data.obj$geno.names[[3]]


	#build a matrix out of the object so we can write it out to a csv file
	final.geno <- rbind(chr, loc, geno.mat)
	final.geno[which(is.na(final.geno))] <- na
	colnames(final.geno) <- marker.names
	
	pheno[which(is.na(pheno))] <- na
	pheno.padding <- matrix(NA, nrow = 2, ncol = dim(pheno)[2])
	final.pheno <- rbind(pheno.padding, pheno)
	
	
	final.data <- cbind(final.pheno, final.geno)
	write.table(final.data, file = filename, sep = ",", quote = FALSE, row.names = FALSE, na = as.character(na))

	
	
}
