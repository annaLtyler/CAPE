#This function flattens the genotype array to a 2D matrix
#for ridge regression, lasso, and elastic net

flatten.geno <- function(geno.obj, ref.allele = "B", verbose = FALSE){
	
	num.total.alleles <- (dim(geno.obj)[2]*dim(geno.obj)[3]) - dim(geno.obj)[3] #subtract out one space from each marker for the reference allele
	
	flat.geno <- matrix(NA, nrow = dim(geno.obj)[1], ncol = num.total.alleles)
	rownames(flat.geno) <- dimnames(geno.obj)[[1]]
	new.colnames <- rep(NA, dim(flat.geno)[2])
	
	allele.names <- dimnames(geno.obj)[[2]]
	ref.locale <- which(allele.names == ref.allele)
	
	start.pos <- 1
	num.alleles <- dim(geno.obj)[2]-1
	
	for(i in 1:dim(geno.obj)[3]){
		if(verbose){report.progress(i, dim(geno.obj)[3])}
		allele.geno <- geno.obj[,-ref.locale,i]
		new.allele.names <- paste(dimnames(geno.obj)[[3]][i], colnames(allele.geno), sep = "_")
		flat.geno[,start.pos:(start.pos+num.alleles-1)] <- allele.geno
		new.colnames[start.pos:(start.pos+num.alleles-1)] <- new.allele.names
		start.pos <- start.pos + num.alleles
		}
	colnames(flat.geno) <- new.colnames
	
	return(flat.geno)
	
	}