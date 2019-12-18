#This script performs the 1D scan of the data
#It first calls the permutation script genome.wide.significance.R
#to calculate the significance threshold. It then uses
#this threshold to determine which markers will be used as
#covariates for each phenotype
#The script creates a table of effect sizes for each marker
#for each phenotype, and a table of covariate flags indicating
#which markers will be used as covariates for each phenotoype
#These are both added to the data object
#The user has the choice to scan the eigentraits (default)
#or the original phenotypes. 
#if the data.object does not specify the genotypes wanted
#in the kinship calculation, the user can add an additional
#object that forces the kinship correction to use more of the
#genome than is represented in the data.obj or geno.obj

singlescan.kin <- function(data.obj, geno.obj = NULL, ref.allele = "A", n.perm = NULL, alpha = c(0.01, 0.05), scan.what = c("eigentraits", "raw.traits"), kin.mats, run.parallel = TRUE, verbose = FALSE, overwrite.alert = TRUE, n.cores = NULL) {

	if(overwrite.alert){
		choice <- readline(prompt = "Please make sure you assign the output of this function to a singlescan.obj, and NOT the data.obj. It will overwrite the data.obj.\nDo you want to continue (y/n) ")
		if(choice == "n"){stop()}
		}

	if(!is.null(kin.mats)){
		use.kinship <- TRUE
		}else{
		use.kinship <- FALSE	
		}
		
	gene <- get.geno(data.obj, geno.obj)
	pheno <- get.pheno(data.obj, scan.what)			
	n.phe = dim(pheno)[2]
	chr.which <- unique(data.obj$chromosome)
	ref.col <- which(colnames(gene) == ref.allele)
	
	#get the covariates and assign the variables
	#to the local environment
	covar.info <- get.covar(data.obj)
	for(i in 1:length(covar.info)){
		assign(names(covar.info)[i], covar.info[[i]])
		}
	covar <- covar.info$covar.names
			

		
	if(is.null(covar)){
		covar.names = "none"
		covar.loc = NULL
		covar.table = NULL
		}

	
	if(is.null(n.perm) || n.perm < 2){alpha = "none"}
	singlescan.obj <- vector(mode = "list", length = 4)
	names(singlescan.obj) <- c("alpha", "alpha.thresh", "covar", "singlescan.results")
	singlescan.obj$covar <- colnames(covar.table)
	singlescan.obj$alpha <- alpha

	#==================================================================
	#if we are using a kinship correction, make sure the phenotypes
	#are mean-centered and there are no missing values in the genotype
	#matrix.
	#==================================================================
	if(use.kinship){
		pheno.means <- apply(pheno, 2, mean)
		tol = 0.01
		non.zero <- intersect(which(pheno.means > 0+tol), which(pheno.means < 0-tol))
		if(length(non.zero) > 0){
			warning("Phenotypes must be mean-centered before performing kinship corrections.")
			cat("Mean-centering phenotypes using norm.pheno()")
			data.obj <- norm.pheno(data.obj)
			}
		missing.vals <- which(is.na(gene))
			if(length(missing.vals) > 0){
				stop("There are missing values in the genotype matrix. Please use impute.geno().")
				}
			}
	#==================================================================


	#==================================================================
	#first do the permutations to get the significance threshold
	#results will be compared to the significance threshold to 
	#determine	
	#==================================================================
	if(!is.null(n.perm) && n.perm >= 2){
		if(verbose){
		cat("\nPerforming permutations to calculate significance threshold...\n")
		}
		# if(run.parallel){
			# alpha.thresh <- genome.wide.threshold.1D.parallel(data.obj, geno.mat = gene, n.perm = n.perm, alpha = alpha, scan.what = scan.what, verbose = verbose, n.cores = n.cores)
			# }else{
		alpha.thresh <- genome.wide.threshold.1D(data.obj, geno.obj = geno.obj, n.perm = n.perm, alpha = alpha, scan.what = scan.what, run.parallel = run.parallel, n.cores = n.cores)
			# }
	singlescan.obj$alpha.thresh <- alpha.thresh
		}else{
		if(verbose){message("Permutations are not being calculated\n\tTo calculate permutations n.perm must be greater than 2.\n")}	
		singlescan.obj$alpha.thresh = "no permutations"
		}
	#==================================================================
	

		#===========================================================================
		#internal functions
		#===========================================================================		

		# #This function gets regression statistics with a
		# #covariate table
		# get.stats <- function(phenoV, markerV, covarV = NULL){
			
			# marker.var <- var(markerV, na.rm = TRUE)
			# if(marker.var == 0){
				# return(rep(NA, 4))
				# }

			# if(use.kinship){
				# if(!is.null(covarV)){
					# model <- lm(phenoV~0+cbind(covarV, markerV))#remove intercept from model	
					# }else{
					# model <- lm(phenoV~0+markerV)#remove intercept from model		
					# }
				# }else{
				# if(!is.null(covarV)){
					# model <- lm(phenoV~cbind(covarV, markerV))
					# }else{
					# model <- lm(phenoV~markerV)
					# }
				# }
			# #take the last line of coefficients.
			# model.coef <- summary(model)$coefficients
			# slope <- model.coef[dim(model.coef)[1],1]
			# se <- model.coef[dim(model.coef)[1],2]
			# t.stat <- abs(model.coef[dim(model.coef)[1],3])
			# p.val <- model.coef[dim(model.coef)[1],4]
							
			# #put together all the statistics we want to keep
			# #we keep the absolute value of the t statistic,
			# #the p value, and the covariate flag
			
			# table.row <- c(slope, se, t.stat, p.val)
			# rm("model.coef", "slope", "se", "t.stat", "p.val", "model")
			# gc()
			# return(table.row)
			# }		
		
		get.chr.stats <- function(ch, phenoV, covar.mat){
			if(verbose){cat("\tCh", ch, "...", sep = "")}
			marker.locale <- which(data.obj$chromosome == ch)
			if(use.kinship){
				
				sink("regress.warnings") #sink regress warnings to a temporary file
				kin.obj <- kinship.on.the.fly(kin.obj = kin.mats, geno = gene, chr1 = ch, chr2 = ch, phenoV = phenoV, covarV = covar.mat)
				sink(NULL) #stop sinking ouput
				
				cor.geno <- kin.obj$corrected.geno
				cor.pheno <- kin.obj$corrected.pheno
				cor.covar <- kin.obj$corrected.covar
				}else{
				cor.geno <- gene
				cor.pheno <- phenoV
				cor.covar <- covar.mat
				}
				
				#calculate the singlescan for the chr markers using the kinship object
				if(run.parallel){
					just.markers <- cor.geno[,marker.locale,drop=FALSE]
					registerDoParallel(cores = n.cores)
					marker.stats <- foreach(m = just.markers, .combine = "cbind") %dopar% {
						get.stats(phenoV = cor.pheno, markerV = m, covarV = cor.covar)
						}
					}else{
	
				marker.stats <- matrix(NA, ncol = length(marker.locale), nrow = 4)
				for(m in 1:length(marker.locale)){
					report.progress(m, length(marker.locale))
					marker.stats[,m]
					test <- get.stats.multiallele(phenotype = cor.pheno, genotype = cor.geno[,,marker.locale[m]], covar.table = cor.covar, ph.family = ph.family, ref.col = ref.col)
					}
				if(verbose){cat("\n")}
		
					}
				if(length(marker.locale) == 1){
					marker.stats <- matrix(marker.stats, ncol = 1)
					}
				colnames(marker.stats) <- data.obj$marker.num[marker.locale]
				return(marker.stats)
				}
		
		
		#===========================================================================
		
		pheno.stats <- vector(mode = "list", length = n.phe)
	
		for(ph in 1:n.phe){
			if(verbose){cat("\nScanning", colnames(pheno)[ph], "\n")}
			phenoV <- pheno[,ph]
			covar.mat <- covar.table
			pheno.stats[[ph]] <- lapply(chr.which, function(x) get.chr.stats(x, phenoV, covar.mat))
			}
	
			
		#reorganize and label pheno.stats
		org.pheno.stats <- vector(mode = "list", length = n.phe)
		names(org.pheno.stats) <- colnames(pheno)
		for(ph in 1:length(pheno.stats)){
			num.markers <- sum(unlist(lapply(pheno.stats[[ph]], function(x) dim(x)[2])))
			stats.mat <- matrix(NA, ncol = 4, nrow = num.markers)
			colnames(stats.mat) <- c("slope", "se", "t.stat", "p.val")
			rownames(stats.mat) <- unlist(lapply(pheno.stats[[ph]], function(x) dimnames(x)[[2]]))
			start.pos <- 1
			for(m in 1:length(pheno.stats[[ph]])){
				nrow.chunk <- dim(pheno.stats[[ph]][[m]])[2]
				stats.mat[start.pos:(start.pos+nrow.chunk-1),] <- t(pheno.stats[[ph]][[m]])
				start.pos <- start.pos+nrow.chunk
				}
			org.pheno.stats[[ph]] <- stats.mat
			}

		pheno.stats <- org.pheno.stats
	

		#calculate statistics for the covariates
		if(!is.null(covar)){
			if(verbose){cat("Calculating statistics for covariates.\n")}
			covar.stats <- vector(mode = "list", length = n.phe)
			for(ph in 1:n.phe){ 
				if(use.kinship){
					sink("regress.warnings") #sink regress warnings to a temporary file
					#calculate the kinship object using the entire genome
					kin.obj <- kinship.on.the.fly(kin.obj = kin.mats, geno = kin.gene.calc, chr1 = NULL, chr2 = NULL, phenoV = pheno[,ph,drop=FALSE], covarV = covar.table)
					sink(NULL) #stop sinking output
					
					cor.geno <- kin.obj$corrected.geno
					cor.pheno <- kin.obj$corrected.pheno
					cor.covar <- kin.obj$corrected.covar
					}else{
					cor.geno <- gene
					cor.pheno <- pheno[,ph]
					cor.covar <- covar.table
					}
					#calculate the singlescan using the covariates
					covar.stats[[ph]] <- apply(matrix(1:dim(cor.covar)[2], nrow = 1), 2, function(x) get.stats(cor.pheno, cor.covar[,x], NULL))
					}

			for(ph in 1:length(pheno.stats)){
				pheno.stats[[ph]] <- rbind(pheno.stats[[ph]], matrix(covar.stats[[ph]], nrow = dim(covar.table)[2], byrow = TRUE))
				rownames(pheno.stats[[ph]])[which(rownames(pheno.stats[[ph]]) == "")] <- colnames(covar.table)
				}
			} #end case for when we have covariates


		unlink("regress.warnings") #delete the temporary file
	
		singlescan.obj$singlescan.results <- pheno.stats
		
		return(singlescan.obj)
	
	}

