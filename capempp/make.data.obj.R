#This function joins phenotype and genoytpe information
#into a cross object. The genotype object going into 
#this function should be retained as the geno.obj
#going forward. This function only adds information
#to the pheno.obj to make it a data.obj
#It also makes sure that the individuals in the
#phenotype and genotypes objects are in the same order

make.data.obj <- function(pheno.obj, geno.obj){
	
	data.obj <- list("pheno" = pheno.obj)
	
	data.obj$marker.num <- 1:length(geno.obj$marker.names)
	data.obj$chromosome <- geno.obj$chromosome
	data.obj$marker.location <- geno.obj$marker.location

	pheno.ind <- rownames(data.obj$pheno)
	geno.ind <- rownames(geno.obj$geno)

	common.ind <- intersect(pheno.ind, geno.ind)
	
	pheno.locale <- match(common.ind, pheno.ind)
	geno.locale <- match(common.ind, geno.ind)

	geno.names <- dimnames(geno.obj$geno)
	names(geno.names) <- c("mouse", "allele", "locus")
	
	# identical(pheno.ind[pheno.locale], geno.ind[geno.locale])
	data.obj$pheno <- pheno.obj[pheno.locale,]
	data.obj$geno.names <- geno.names
	geno.obj$geno <- geno.obj$geno[geno.locale,,]
	
	
	return(list("data.obj" = data.obj, "geno.obj" = geno.obj$geno))

	
}