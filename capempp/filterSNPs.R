#This function uses matrixEpistasis to prioritize SNPs for 
#testing in cape
#covariates that are included in the data.obj are automatically
#included here.
#the tests are performed in chunks of 10k snps
# data.obj <- readRDS("~/Documents/Data/Scleroderma/project_data/cross.RData")
# geno.obj <- readRDS('~/Documents/Data/Scleroderma/project_data/geno.RData')


filterSNPs <- function(data.obj = NULL, geno.obj = NULL, target.markers = 1500, ref.allele = "A", pval = 0.01, chunk.size = 5000, run.matrixEpistasis = TRUE, n.cores = 4, run.parallel = TRUE){

	snp.files <- list.files(pattern = "SNP.Filter.Results")
	if(length(snp.files) == 0){run.matrixEpistasis <- TRUE}

	if(run.matrixEpistasis){
		require(MatrixEpistasis)


		pheno <- data.obj$pheno
		geno <- get.geno(data.obj, geno.obj)

		#get the covariates and assign the variables
		#to the local environment
		covar.info <- get.covar(data.obj)
		for(i in 1:length(covar.info)){
			assign(names(covar.info)[i], covar.info[[i]])
			}

		ref.locale <- which(dimnames(geno)[[2]] == ref.allele)
		the.rest <- setdiff(1:dim(geno)[[2]], ref.locale)
		geno.mat <- geno[,the.rest,]

		chunk.geno <- chunkV(1:ncol(geno.mat), ceiling(ncol(geno.mat)/chunk.size))
		chunk.pairs <- pair.matrix(1:length(chunk.geno), self.pairs = TRUE)
		}

	#====================================================
	#internal functions
	#====================================================
	get.sig.snps <- function(pair.ind, pheno.mat){
		snp1.ind <- chunk.geno[[chunk.pairs[pair.ind,1]]]
		snp2.ind <- chunk.geno[[chunk.pairs[pair.ind,2]]]
		snp1 <- geno.mat[,snp1.ind]
		snp2 <- geno.mat[,snp2.ind]
	
		res <- matrixEpistasis(snp1, snp2, pheno.mat, covariate = covar.table)
		r = res$r 
		df = res$df	
		
		corrThreshold <- p2c(pval = pval, df)
		index <- which(abs(r) > corrThreshold, arr.ind=TRUE )
	
		# get the SNP names
		sig.snp1 <- colnames(snp1)[index[,1]]
		sig.snp2 <- colnames(snp2)[index[,2]]
	
		# use matrixPval function to calculate p values for only those of interest
		pvalue <- matrixPval(r[index] ,df)
	
		sig.table <- data.frame(sig.snp1, sig.snp2, pvalue)
		return(sig.table)
		
		}	

	merge.pairs <- function(tableX){
		sorted.pairs <- t(apply(tableX, 1, function(x) sort(x[1:2])))
		pair.text <- apply(sorted.pairs, 1, function(x) paste(x[1], x[2], sep = "_"))
		return(pair.text)
		}

	pheno.which <- function(marker){
		test.pheno <- lapply(snp.lists, function(x) which(x[,1:2] == marker))
		is.present <- unlist(lapply(test.pheno, length))
		return(is.present)
		}
	#====================================================

		
	if(run.matrixEpistasis){

		cat("\n")
		for(i in 1:ncol(pheno)){
			cat("Analyzing", colnames(pheno)[i], "...\n")

			if (run.parallel) {
				cl <- makeCluster(n.cores)
				registerDoParallel(cl)
				pheno.sig <- foreach(p = 1:nrow(chunk.pairs), .combine = "rbind", .export = c("matrixEpistasis", "p2c", "matrixPval")) %dopar% {
					get.sig.snps(p, pheno[,i,drop=FALSE])
				}
				stopCluster(cl)
			} else {

				pheno.sig <- c()
				index <- 1:nrow(chunk.pairs)
				for (p in index) {
				  pheno.sig <- rbind(pheno.sig, get.sig.snps(p, pheno[,i,drop=FALSE]))
				}

		  }
			saveRDS(pheno.sig, paste0("SNP.Filter.Results.", colnames(pheno)[i], ".RData"))
			cat(nrow(pheno.sig), '\tsignificant pairs found\n')
		}
		snp.files <- list.files(pattern = "SNP.Filter.Results")
	}
	
	file.names <- multi.strsplit(snp.files, patterns = c("SNP.Filter.Results.", ".RData"))
	
	snp.lists <- lapply(snp.files, readRDS)

	#Now we will go through some machinations to get SNPs that
	#influence as many traits as possible
	
	#Find pairs of SNPs that influence pairs of traits
	merged.snp.pairs <- lapply(snp.lists, merge.pairs)
	plotVenn(merged.snp.pairs, cat.names = file.names, lwd = 3)

	#now find all pairs that influence at least pairs of traits.
	pheno.pairs <- pair.matrix(1:length(snp.lists))
	all.shared <- apply(pheno.pairs, 1, function(x) intersect(merged.snp.pairs[[x[1]]], merged.snp.pairs[[x[2]]]))

	total.snps <- unique(unlist(lapply(unlist(all.shared), function(x) strsplit(x, "_"))))

	#if we don't have enough SNPs yet start looking for SNPs
	#that influence multiple traits
	if(length(total.snps) < target.markers){
		#get the names of all SNPs influencing all phenotypes
		ind.snps <- unique(unlist(lapply(snp.lists, function(x) unique(c(unique(as.character(x[,1])), unique(as.character(x[,2])))))))

		if(!file.exists("SNP.pairs.in.Pheno.RData")){
			snp.pheno <- Reduce("rbind", lapply_pb(ind.snps, pheno.which))
			saveRDS(snp.pheno, "SNP.pairs.in.Pheno.RData")
			}else{
			snp.pheno <- readRDS("SNP.pairs.in.Pheno.RData")
			}

		num.pheno <- apply(snp.pheno, 1, function(x) length(which(x > 0)))

		num.affected <- length(snp.lists)
		while(length(total.snps) < target.markers && num.affected > 0){
			pheno.affected.locale <- which(num.pheno == num.affected)
			snp.traits <- ind.snps[pheno.affected.locale]
			test.list <- unique(c(total.snps, snp.traits))

			#if all SNPs will fit in the list, just add them
			if(length(test.list) <= target.markers){
				total.snps <- test.list
				num.affected <- num.affected - 1
				}else{
				#otherwise, add to the list based on how many influences each SNP has
				sub.pheno <- snp.pheno[pheno.affected.locale,]
				sub.snps <- ind.snps[pheno.affected.locale]
				total.influence <- rowSums(sub.pheno)
				sub.pheno <- sub.pheno[order(total.influence, decreasing = TRUE),]
				sub.snps <- sub.snps[order(total.influence, decreasing = TRUE)]
				num.to.take <- target.markers - length(total.snps)
				while(length(total.snps) < target.markers && length(sub.snps) > 0){
					total.snps <- unique(c(total.snps, sub.snps[1:num.to.take]))
					sub.pheno <- sub.pheno[-(1:num.to.take),,drop=FALSE]
					sub.snps <- sub.snps[-(1:num.to.take)]
					num.to.take <- target.markers - length(total.snps)
					}
				}
			}
		}

	write.table(total.snps, "filteredSNPs.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
		
	invisible(total.snps)
		
	}
