#This function selects phenotypes in a cross
#object. 
#min.entries is the minimum number of data entries
#the phenotype needs to have for it to be valid
#If scale.pheno is set to TRUE, the phenotypes are 
#mean centered and standarized.
#if rank.norm.pheno is TRUE, the phenotypes are 
#rank Z normalized
#only numeric phenotypes can be used.

select.pheno <- function(data.obj, pheno.which, min.entries = 5, scale.pheno = FALSE, rank.norm.pheno = FALSE){
	check.underscore(data.obj)
	# check.bad.markers(data.obj)
	
	
	pheno <- data.obj$pheno

	#find the phenotype column numbers if 
	#names have been put in instead of numbers	
	pheno.num <- get.col.num(pheno, pheno.which)
	
	new.pheno <- pheno[,pheno.num]

		# #make sure the phenotypes are numeric
		# #and replace the phenotype matrix with 
		# #the selected phenotypes
		new.pheno <- apply(new.pheno, 2, as.numeric)
		
		#check to see if there are any phenotypes with
		#fewer than 5 entries
		data.entries <- as.vector(apply(new.pheno, 2, function(x) length(which(!is.na(x)))))
		bad.pheno <- which(data.entries <= min.entries)
		
		if(length(bad.pheno) > 0){
			final.pheno <- new.pheno[,-bad.pheno]
			message("The following phenotypes had fewer than ", min.entries, " entries and were removed.\n")
			cat(paste("\t", colnames(new.pheno)[bad.pheno], "\n"))
			}else{
				final.pheno <- new.pheno
				}
				
	#This function mean centers and standardizes a vector
	center.std <- function(v){
		mean.v <- mean(v, na.rm = TRUE)
		centered <- v - mean.v
		sd.v <- sd(v, na.rm = TRUE)
		final.v <- centered/sd.v
		return(final.v)
		}

	if(rank.norm.pheno){
		final.pheno <- apply(final.pheno, 2, rz.transform)
		}

	if(scale.pheno){
		final.pheno <- apply(final.pheno, 2, center.std) #mean center and standardize the phenotypes
		}

		rownames(final.pheno) <- rownames(pheno)		
		data.obj$pheno <- final.pheno
			
	return(data.obj)
	
	}