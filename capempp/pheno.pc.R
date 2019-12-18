#This function breaks phenotypes into subsets
#and creates a new set of phenotypes based on
#the first principal component of each subgroup
#make sure there are no missing data before
#running this 
pheno.pc <- function(data.obj, pheno.list){
	
	
	pheno <- data.obj$pheno

	if(length(which(is.na(pheno))) > 0){
		stop("There cannot be missing data in the phenotype matrix.")
		}
	
	pc.pheno <- matrix(NA, nrow = nrow(pheno), ncol = length(pheno.list))
	colnames(pc.pheno) <- names(pheno.list)
	for(i in 1:length(pheno.list)){
		pheno.locale <- which(colnames(pheno) %in% pheno.list[[i]])
		sub.pheno <- pheno[,pheno.locale]
		pcs <- svd(sub.pheno)
		pc.pheno[,i] <- pcs$u[,1]	
		}
	return(pc.pheno)
	
}