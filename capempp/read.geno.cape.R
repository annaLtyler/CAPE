#This function reads in a data file in the r/qtl format
#It converts letter genotypes to numbers if required.
#It parses the data into a genotype object
#if filename is left empty, the script will ask the
#use to choose a file.
#This function is useful when the genotype matrix 
#is too large to reasonably contain in the data.obj
#after using read.geno and read.pheno, make.data.obj
#should be used to combine information from the 
#genotype and phenotype objects into a data.obj. 
#After this step the result of make.data.obj is
#the data.obj, and the genotype object remains
#unchanged as the genotype object. Both need to
#be used in subsequent functions.
#The phenotype object can be discarded or overwritten.


read.geno.cape <- function(file.format = c("cape", "csv", "rdata"), filename = NULL, marker.info.file = NULL, geno.col = NULL, delim = ",", na.strings = "-", check.chr.order = TRUE, allele.names = NULL) {
	require(abind)

	if(file.format == "csv"){
		if(is.null(filename)){
			filename <- file.choose()
			}
			
		cross.data <- read.table(filename, na.strings = na.strings, stringsAsFactors = FALSE, sep = delim, header = TRUE)
	
	
	    #determine where phenotypes end and genotypes begin by blanks in row 1
		beginGeno = match(FALSE, is.na(suppressWarnings(as.numeric(cross.data[1,]))))
	
		chr <- as.vector(as.matrix(cross.data[1,beginGeno:dim(cross.data)[2]]))
		
		if(check.chr.order){
			x.locale <- grep("x", chr, ignore.case = TRUE)
			y.locale <- grep("y", chr, ignore.case = TRUE)
			just.num.chr <- setdiff(1:length(chr), c(x.locale, y.locale))
			consec.chr <- consec.pairs(as.numeric(chr[just.num.chr]))
			order.check <- apply(consec.chr, 1, function(x) x[2] - x[1])
			if(length(which(order.check < 0)) > 0){
				warning("The chromosomes appear to be out of order.\nIt is best to sort the chromosomes before beginning the analysis.")
				}
			}
			
		marker.loc <- as.numeric(cross.data[2,beginGeno:dim(cross.data)[2]])
	
		#take out the genotype matrix
		#It begins in the 3rd row after chromosome numbers and marker locations
		geno <- as.matrix(cross.data[3:dim(cross.data)[1],beginGeno:dim(cross.data)[2]])
	   
	   if(is.null(geno.col)){
		   	geno.col <- 1:dim(geno)[2]
	   		}
	   	geno.columns <- get.col.num(geno, geno.col)
	  	geno <- geno[,geno.columns] 
	  	chr <- chr[geno.columns]
		marker.loc <- marker.loc[geno.columns]
	   
	    #run a check to see how the genotypes are stored.
		#genotypes can be stored as (0,1,2), ("A","H","B")
		#or as probabilities between 0 and 1
		#if the genotypes are stored as (0,1,2) or 
		#letters, we need to convert them to probabilities
	    
	    found.genotype.mode <- 0 #create a flag to determine whether we have figured out how the genotypes are encoded
	    
	    genotype.class <- class(cross.data[,beginGeno])
	    if(genotype.class != "character"){
		    all.genotypes <- sort(unique(na.omit(as.numeric(as.matrix(geno))))) #get a vector of all genotypes used
		    }else{
			    all.genotypes <- sort(unique(na.omit(as.vector(as.matrix(geno)))))
			    if(length(all.genotypes) > 3){
			    	#look for empty genotypes
			    	cat("I have detected", length(all.genotypes), "genotypes:\n")
			    	print(all.genotypes)
					stop("Please check for missing genotype values or other errors in the genotype data.")
			    	}
		    	}
		
		#check to see if the genotypes are encoded as letters
	
	  	if(genotype.class == "character"){
	  		het.present <- grep("H", all.genotypes)
	  		if(length(all.genotypes) > 2 && length(het.present) == 0){
	  			stop("Heterozygotes must be coded by H")
	  			}
	  			
	 		
	  		found.genotype.mode <- 1 #indicate that we have found the genotype mode for the file
			#assign baseGeno and notBaseGeno
		    #the baseGeno is assigned a numeric value
		    #of 0. The notBaseGeno is assigned 1
		    #by default we make the first letter
		    #alphabetically the base genotype
		    if(length(all.genotypes) == 3){
				baseGeno <- all.genotypes[all.genotypes != "H"][1]
				notBaseGeno <- all.genotypes[all.genotypes != "H"][2]
				cat("The genotypes are encoded as ", baseGeno, ", H, ", notBaseGeno, "\nConverting to 0, 0.5, 1.\n", sep = "")
				}else{
				# baseGeno <- all.genotypes[all.genotypes != "H"][1]
				# notBaseGeno <- "H"
				baseGeno <- sort(all.genotypes)[1]
				notBaseGeno <- sort(all.genotypes)[2]
				cat("The genotypes are encoded as ", baseGeno, " and ", notBaseGeno, "\nConverting to 0 and 1.\n", sep = "")
		 		}
		
	
		
		    #turn baseGeno, H, and notBaseGeno to 0, .5, and 1 respectively
			#This function takes in a vector and converts the letters to
			#the appropriate numbers
		
			convert.geno.letter <- function(genotypes){
				genotypes[which(as.character(genotypes) == baseGeno)] <- 0
		    	genotypes[which(as.character(genotypes) == notBaseGeno)] <- 1
		    	if(length(all.genotypes) == 3){
		    		genotypes[which(as.character(genotypes) == "H")] <- 0.5
		    		}else{
		    		genotypes[which(as.character(genotypes) == "H")] <- 1
		    		}
		    	return(as.numeric(genotypes))
		    	}
		 	geno <- apply(geno, 2, convert.geno.letter) 
	  		}
	    
		
	 	#check to see if the genotypes are encoded as (0, 1, 2)
	 	numeric.test <- which(all.genotypes == 2) #check for 2, since 2 is unique to this encoding
	 	if(length(numeric.test) > 0){
	 		outside.upper.bound <- which(all.genotypes > 2)
	 		outside.lower.bound <- which(all.genotypes < 0)
	 		if(length(outside.upper.bound) > 0 || length(outside.lower.bound) > 0){
	 			stop("Assuming (0,1,2) coding, but I detected genotypes greater than 2 or less than 0.")
	 			}
	 		cat("The genotypes are encoded as 0, 1, 2.\nConverting to 0, 0.5, 1.\n")
	 		found.genotype.mode <- 1 #set the flag indicating we've figured out the encoding
	 		#turn 0, 1, 2 into 0, 0.5 and 1 respectively
			convert.geno.number <- function(genotypes){
		        genotypes[which(as.numeric(genotypes) == 1)] <- 0.5
		        genotypes[which(as.numeric(genotypes) == 2)] <- 1
		        return(as.numeric(genotypes))
				}
		
		 	geno <- apply(geno, 2, convert.geno.number) 
	
		 	}
	 	
	 	#if we still haven't found the genotype mode yet
	 	#check to see if the genotypes are encoded as probabilities
	 	if(found.genotype.mode == 0){
	 		min.geno <- min(all.genotypes); max.geno <- max(all.genotypes) #find the max and min values for genotypes
	 		if(min.geno >= 0 && max.geno <= 1){ #and make sure they are bounded as probabilities are
	 			found.genotype.mode <- 1 #set the flag to indicate we have found the genotype mode
	 			cat("The genotypes are encoded as probabilities.\nNo conversion needed.\n")
	 			
				#we still need to conver the data frame to a numeric matrix for later compatibility
	 			convert.geno.prob <- function(genotypes){
			        return(as.numeric(genotypes))
					}
		
			 	geno <- apply(geno, 2, convert.geno.prob) 
	
	 			}
	 		}
	 	
	 	
	 	#If after all this, we haven't found the genotype encoding, stop and warn the user
	 	if(found.genotype.mode == 0){
	 		stop("\nGenotypes must be encoded as (0, 1, 2), (A,H,B), or probabilities.\n")
	 		}
	
	
	
		#take out the sex chromosomes and invariant markers
		x.locale <- grep("X", chr, ignore.case = TRUE)
		if(length(x.locale) > 0){
			message("\nRemoving markers on the X chromosome")
			geno <- geno[,-x.locale]
			chr <- chr[-x.locale]
			marker.loc <- marker.loc[-x.locale]
			}
			
		y.locale <- grep("Y", chr, ignore.case = TRUE)
		if(length(y.locale) > 0){
			message("\nRemoving markers on the Y chromosome")
			geno <- geno[,-y.locale]
			chr <- chr[-y.locale]
			marker.loc <- marker.loc[-y.locale]
			}
			
		#take out markers with only 1 allele
		num.allele <- apply(geno, 2, function(x) length(unique(x)))
		mono.allele <- which(num.allele == 1)
		if(length(mono.allele) > 0){
			message("\nRemoving invariant markers.\n")
			geno <- geno[,-mono.allele]
			chr <- chr[-mono.allele]
			marker.loc <- marker.loc[-mono.allele]
			}
		
	
	
		na.locale <- which(is.na(geno))
		if(length(na.locale) > 0){
			message("I have detected missing values in the genotype matrix.\n\tIf you are planning to use the kinship correction, please use impute.geno() to impute the genotype data.\n")
			}
	
		#put in code here to distribute the genotypes between -1 and 1 so we get symmetric m12/m21 null distributions
		#construct the data object
		marker.names <- colnames(geno)
		colnames(geno) <- 1:dim(geno)[2]
		rownames(geno) <- 1:dim(geno)[1]
	
		#scale the genotypes to lie between 0 and 1
		#even if it's a backcross
		geno <- geno/max(geno, na.rm = TRUE)
		marker.num <- 1:dim(geno)[2]
	
		final.data <- list(geno, chr, marker.names, marker.num, marker.loc)
		names(final.data) <- c("geno", "chromosome", "marker.names", "marker.num", "marker.location")
		     
		cat("Read in the following data:\n")
		cat("\t-", dim(geno)[1], "individuals -\n")
		cat("\t-", dim(geno)[2], "markers -\n")
	
	    return(final.data)
		}    
	
	
		if(file.format == "rdata"){
			if(is.null(filename)){
				cat("Please choose the genotype file.\n")
				filename <- file.choose()
				}
		geno <- readRDS(filename)
		
		if(is.null(marker.info.file)){
			cat("Please choose the marker info file.\n")
			marker.info.file <- file.choose()
			}

		marker.info <- read.table(marker.info.file, sep = delim, header = TRUE, stringsAsFactors = FALSE)
	
		#make sure marker names in the info file and genotype object line up
		marker.dim = length(dim(geno))
		
		names.same <- identical(marker.info[,1], dimnames(geno)[[marker.dim]])	
		if(!names.same){
			warning("Marker names in the genotype object and info file are not identical. \nFiltering to only common markers...\n")
			common.markers <- intersect(marker.info[,1], dimnames(geno)[[marker.dim]])
			info.locale <- match(common.markers, marker.info[,1])
			geno.locale <- match(common.markers, dimnames(geno)[[marker.dim]])
			marker.info <- marker.info[info.locale,]
			if(marker.dim == 3){
				geno <- geno[,,geno.locale]
				}else{
				geno <- geno[,geno.locale]	
				}
			}
		
		#take out the x and y chromosomes (and chromosome 0), if these
		#are present in marker.info since we will not test these
		sex.locale <- c(which(marker.info[,2] == "X"), which(marker.info[,2] == "Y"), which(marker.info[,2] == 0))
		if(length(sex.locale) > 0){
			if(marker.dim == 3){
				geno <- geno[,,-sex.locale]
			}else{
				geno <- geno[,-sex.locale]
				}
			}
			
		#sort the markers by chromosome then position
		sorted.snps <- sortByThenBy(marker.info, sort.cols = c(2,3), col.type = c("n", "n"))
	
		#make sure the genotypes are in the same order
		if(!identical(marker.info[,1], sorted.snps[,1])){
			snp.order <- match(sorted.snps[,1], dimnames(geno)[[marker.dim]])
			if(marker.dim == 3){
				geno <- geno[,,snp.order]
				}else{
				geno <- geno[,snp.order]	
				}
			}
			
		#finally turn a two-dimensional matrix into a 3D array
		if(marker.dim == 2){
			misordered <- abind(geno, 1-geno, along = 3)
			geno <- aperm(misordered, c(1,3,2))
			if(is.null(allele.names)){
				dimnames(geno)[[2]] <- c("A", "B")
				}else{
				dimnames(geno)[[2]] <- allele.names	
				}
			}
					
		results <- list(geno, sorted.snps[,1], sorted.snps[,2], sorted.snps[,3], 1:nrow(sorted.snps))
		names(results) <- c("geno", "marker.names", "chromosome", "marker.location", "marker.num")	
		
		return(results)
		
	    }
}
