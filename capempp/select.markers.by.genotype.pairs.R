
select.markers.by.genotype.pairs <- function(data.obj, geno.obj, ref.allele = "B", max.markers.per.chunk = 5000, min.per.genotype = 10, n.cores = 4, run.parallel = TRUE){

	
	geno <- get.geno(data.obj, geno.obj)
	ref.locale <- which(colnames(geno) == ref.allele)
	if(length(ref.locale) == 0){
		stop("can't find reference allele.")
		}
	geno <- geno[,-ref.locale,]
	
	#=======================================================
	#chunk up the genome so our vector and matrices don't 
	#get too big
	#=======================================================
	
	num.markers <- dim(geno)[[3]]
	num.bins <- ceiling(num.markers/max.markers.per.chunk)
	chunked.markers <- chunkV(1:dim(geno)[[3]], num.bins)
	chunk.pairs <- pair.matrix(1:length(chunked.markers), self.pairs = TRUE)
	
	#=======================================================
	# internal functions
	#=======================================================
	
	bin.full.marker <- function(full.marker){
		binned.marker <- apply(full.marker, 2, bin.vector)
		return(binned.marker)
		}
	
	
	count.geno <- function(marker1.geno, marker2.geno, min.per.genotype){
		min.per.pair <- apply(allele.pairs, 1, function(x) min(table(marker1.geno[,x[1]], marker2.geno[,x[2]])))
		good.locale <- which(min.per.pair > min.per.genotype)
		if(length(good.locale) > 0){
			final.result <- cbind(allele.pairs[good.locale,,drop=FALSE], min.per.pair[good.locale])
			return(final.result)		
			}
		}
		
	count.geno.by.bin <- function(pair.idx){
		geno1.idx <- chunked.markers[[chunk.pairs[pair.idx,1]]]
		geno2.idx <- chunked.markers[[chunk.pairs[pair.idx,2]]]
		all.marker.pairs <- cbind(rep(geno1.idx, length(geno2.idx)), rep(geno2.idx, each = length(geno1.idx)))
		pair.names <- apply(all.marker.pairs, 1, function(x) paste(x[1], x[2], sep = "_"))

		# test.num = 10;warning("I'm running a small number for testing purposes")
		test.num <- nrow(all.marker.pairs)

		pair.results <- lapply(1:test.num, function(x) count.geno(marker1.geno = binned.mat[,,all.marker.pairs[x,1]], marker2.geno = binned.mat[,, all.marker.pairs[x,2]], min.per.genotype = min.per.genotype))
		
		names(pair.results) <- pair.names[1:test.num]
		not.null <- unlist(lapply(pair.results, function(x) !is.null(x)))
		if(length(which(not.null)) > 0){
			not.null.pairs <- pair.results[which(not.null)]
			for(i in 1:length(not.null.pairs)){
				rownames(not.null.pairs[[i]]) <- rep(names(not.null.pairs)[i], nrow(not.null.pairs[[i]]))
				}
			not.null.mat <- Reduce("rbind", not.null.pairs)
			return(not.null.mat)
			}
		}
	#=======================================================
	
	
	#=======================================================
	#bin the genotype matrix into 0, 0.5, and 1
	#and reassemble
	#=======================================================	
	cat("binning genotypes into 0, 0.5, and 1...\n")
	binned.geno <- lapply_pb(1:dim(geno)[[3]], function(x) bin.full.marker(geno[,,x]))
	binned.mat <- abind(binned.geno, along = 3)
	alleles <- colnames(geno)
	allele.pairs <- pair.matrix(alleles, self.pairs = TRUE, ordered = TRUE)
	
	geno.pairs <- pair.matrix(1:dim(geno)[[3]])
	
	#=======================================================
	#count up the animals in all genotype pair bins
	#=======================================================	
	cat("testing genotype pairs...\n")

	if (run.parallel) {

	  cl <- makeCluster(n.cores)
	  registerDoParallel(cl)
	  allele.table <- foreach(p = c(1:nrow(chunk.pairs)), .combine = "rbind") %dopar% {
	    count.geno.by.bin(p)
	  }
	  stopCluster(cl)

	} else {

	  allele.table <- c()
	  index <- 1:nrow(chunk.pairs)
	  for (p in index) {
	    allele.table <- rbind(allele.table, count.geno.by.bin(p))
	  }

	}
	
	split.markers <- strsplit(rownames(allele.table), "_")
	marker1 <- as.numeric(unlist(lapply(split.markers, function(x) x[1])))
	marker2 <- as.numeric(unlist(lapply(split.markers, function(x) x[2])))
	
	allele.table <- cbind(marker1, marker2, allele.table)
	colnames(allele.table) <- c("marker1", "marker2", "allele1", "allele2", "num.ind")
	
	write.table(allele.table, "Allele.Pairs.With.Representation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

	#also rebuild the genotype and data objects to keep only the alleles that are in 
	allele1 <- match(allele.table[,3], colnames(geno))
	allele2 <- match(allele.table[,4], colnames(geno))
	
	u_marker1 <- unique(cbind(marker1, allele1))
	u_marker2 <- unique(cbind(marker2, allele2))
	u_markers <- unique(rbind(u_marker1, u_marker2))
	
	new.geno <- apply(u_markers, 1, function(x) geno[,x[2],x[1]])
	new.geno.names <- apply(u_markers, 1, function(x) paste(dimnames(geno)[[3]][x[1]], dimnames(geno)[[2]][x[2]], sep = "."))
	colnames(new.geno) <- new.geno.names
	
	new.geno.obj <- abind(1-new.geno, new.geno, along = 1.5)
	all.alleles <- c("A", "B")
	not.ref <- setdiff(all.alleles, ref.allele)
	colnames(new.geno.obj) <- c(ref.allele, not.ref[1])
	
	#now we need to update the marker information in the data.obj to 
	#reflect the new genotype object
	data.obj$geno.names$allele <- colnames(new.geno.obj)
	
	just.markers <- dimnames(geno)[[3]][u_markers[,1]]
	new.chr <- get.marker.chr(data.obj, just.markers)
	new.loc <- get.marker.location(data.obj, just.markers)
	new.num <- get.marker.num(data.obj, just.markers)
	
	marker.order <- order(new.num)
	new.geno.obj <- new.geno.obj[,,marker.order]
	
	data.obj$chromosome <- new.chr[marker.order]
	data.obj$marker.location <- new.loc[marker.order]
	data.obj$marker.num <- new.num[marker.order]
	data.obj$geno.names$locus <- dimnames(new.geno.obj)[[3]][marker.order]

	names(dimnames(new.geno.obj)) <- c("mouse", "allele", "locus")

	final.results <- list("data.obj" = data.obj, "geno.obj" = new.geno.obj)
	return(final.results)
	
	}
	
	
	
	
	
	
	
	
