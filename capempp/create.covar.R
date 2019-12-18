#This script takes a variable from the phenotype matrix
#for example, diet treament or sex and creates a marker
#variable that can be used as a covariate.
#It creates a marker that is numeric and assigns the 
#numeric value to each of the allels at all loci for 
#the given individual.


create.covar <- function(data.obj, pheno.which){

	check.underscore(data.obj)
	# check.bad.markers(data.obj)
	
	pheno <- data.obj$pheno
	# geno <- get.geno(data.obj, geno.obj)
	# chromosome <- data.obj$chromosome
	# marker.location <- data.obj$marker.location
	
	marker.locale <- get.col.num(pheno, pheno.which)

	#make a separate covariate table, then modify the dimnames
	#in the genotype object to include the covariates
	#do not modify the genotype object
	
	covar.table <- pheno[,marker.locale,drop=FALSE]
	rownames(covar.table) <- rownames(pheno)
	data.obj$covar.table <- covar.table
	
		
	#take the phenotypes made into markers out of the phenotype matrix
	new.pheno <- pheno[,-marker.locale]
	data.obj$pheno <- new.pheno
	data.obj$non.allelic.covar <- pheno.which
	data.obj$geno.names[[3]] <- c(data.obj$geno.names[[3]], pheno.which)
	data.obj$chromosome <- c(data.obj$chromosome, rep(0, length(pheno.which)))
	data.obj$marker.location <- c(data.obj$marker.location, 1:length(pheno.which))
	
	return(data.obj)
	
	}