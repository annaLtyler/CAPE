#This function uses significant SNPs to define linkage blocks


linkage.blocks <- function(data.obj, geno.obj = NULL, p.or.q = 0.05, collapse.linked.markers = TRUE, r2.thresh = 0.8){
	
	geno.names <- data.obj$geno.names
	marker.names <- geno.names[[3]]
	net.data <- data.obj$var.to.var.p.val
	pheno.net.data <- data.obj$max.var.to.pheno.influence
	data.obj$network.p.or.q <- p.or.q
	
	if(collapse.linked.markers){
		data.obj$r2.thresh <- r2.thresh
		}

	if(length(net.data) == 0){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}
	
	#get the markers with significant influences on other markers
	#and on phenotypes
	var.sig.col <- which(colnames(var.influences) == "p.adjusted")
	sig.markers <- net.data[which(net.data[, var.sig.col] <= p.or.q), 1:2]
	
	sig.pheno <- vector(mode = "list", length = length(pheno.net.data))
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.net.data[[1]]))))
	for(ph in 1:length(pheno.net.data)){
		sig.locale <- which(as.numeric(pheno.net.data[[ph]][,pheno.sig.col]) <= p.or.q)
		sig.pheno[[ph]] <- pheno.net.data[[ph]][sig.locale,1]
		}

	if(length(c(sig.markers, unlist(sig.pheno))) == 0){
		stop("There are no significant markers at p = ", p.or.q)
		}

	#find all the markers with significant p values
	#and sort them
	u_markers <- unique(c(sig.markers, unlist(sig.pheno)))
	just.markers <- sapply(strsplit(u_markers, "_"), function(x) x[[1]][1])
	marker.locale <- match(just.markers, marker.names)
	sorted.markers <- u_markers[order(marker.locale)]
	marker.chr <- data.obj$chromosome[sort(marker.locale)]
	u_chr <- unique(marker.chr)
	
	#take these out of the genotype matrix
	geno <- get.geno(data.obj, geno.obj)
	marker.geno <- geno[,,marker.locale]
	
	get.cor <- function(marker1, marker2){
		marker1.split <- strsplit(marker1, "_"); marker1.name <- marker1.split[[1]][1]; marker1.allele <- marker1.split[[1]][2]
		marker2.split <- strsplit(marker2, "_"); marker2.name <- marker2.split[[1]][1]; marker2.allele <- marker1.split[[1]][2]
		if(is.na(marker1.allele)){marker1.allele <- 1}
		if(is.na(marker2.allele)){marker2.allele <- 1}
		marker.cor <- (cor(marker.geno[,marker1.allele,marker1.name], marker.geno[,marker2.allele,marker2.name], use = "pairwise.complete.obs"))^2

		return(marker.cor)
		}

	#go through each chromosome separately and find the linkage blocks on each chromosome
	linkage.blocks <- list()
	num.blocks <- 1
	for(ch in 1:length(u_chr)){
		chr.markers <- sorted.markers[which(marker.chr == u_chr[ch])]

		if(!collapse.linked.markers){
			for(i in 1:length(chr.markers)){
				linkage.blocks[[num.blocks]] <- chr.markers[i]
				names(linkage.blocks)[num.blocks] <- chr.markers[i]
				num.blocks <- num.blocks + 1
				}
			}else{
		
		
		if(length(chr.markers) == 1){
			chr.allele <- strsplit(chr.markers, "_")[[1]][2]
			linkage.blocks[[num.blocks]] <- chr.markers
			names(linkage.blocks)[num.blocks] <- paste("Chr", u_chr[ch], "_", 1, "_", chr.allele, sep = "")
			num.blocks <- num.blocks + 1
			}else{
		
			# just.markers <- unique(sapply(strsplit(chr.markers, "_"), function(x) x[[1]][1]))
			marker.pairs <- pair.matrix(chr.markers)
			
			all.cor <- matrix(NA, nrow = length(chr.markers), ncol = length(chr.markers))
			rownames(all.cor) <- colnames(all.cor) <- chr.markers
	
			for(i in 1:length(marker.pairs[,1])){
				all.cor[as.character(marker.pairs[i,1]), as.character(marker.pairs[i,2])] <- get.cor(marker1 = marker.pairs[i,1], marker2 = marker.pairs[i,2])
				}

			#zero out the diagonal and the lower triangle
			all.cor[lower.tri(all.cor, diag = TRUE)] <- 0
			in.ld <- which(all.cor >= r2.thresh, arr.ind = TRUE)


			#if there are no markers in LD, add each of the 
			#markers to the total blocks for the whole genome
			if(length(in.ld) == 0){
				for(i in 1:length(chr.markers)){
					linkage.blocks[[num.blocks]] <- chr.markers[[i]]
					names(linkage.blocks)[num.blocks] <- chr.markers[[i]]
					num.blocks <- num.blocks + 1
					}
				}else{
					
				#otherwise, go through the in.ld matrix and combine markers that are linked
				linked.markers <- list(in.ld[1,])
	
				if(dim(in.ld)[1] > 1){
					for(i in 2:length(in.ld[,1])){

						row.to.check <- in.ld[i,]
						#look for blocks that already contain markers we're looking at
						#The lapply function subtracts the length of the combined unique
						#vector from the length of the total vector. If the result is 
						#positive, there are common nodes between the new row and an
						#existing block
						common.nodes <- lapply(linked.markers, function(x) length(c(x,row.to.check))-length(unique(c(x,row.to.check))))
						shared.node.locale <- which(common.nodes > 0)

						#if there are no shared nodes, start a new block 
						if(length(shared.node.locale) == 0){
							linked.markers[[(length(linked.markers)+1)]] <- as.numeric(row.to.check)
							}

						#if we need to combine the new row with one other block
						if(length(shared.node.locale) == 1){
							linked.markers[[shared.node.locale]] <- as.numeric(unique(c(linked.markers[[shared.node.locale]], row.to.check)))
							}
		
						if(length(shared.node.locale) > 1){
							stop("There is more than one common node")
							}

						}
					}

					#add these blocks to the total list
					for(i in 1:length(linked.markers)){
						marker.names <- chr.markers[linked.markers[[i]]]
						marker.alleles <- sort(unique(unlist(lapply(strsplit(chr.markers[linked.markers[[i]]], "_"), function(x) x[2]))))
						marker.allele.text <- paste(marker.alleles, collapse = "_")
						linkage.blocks[[num.blocks]] <- marker.names
						names(linkage.blocks)[num.blocks] <- paste("Chr", u_chr[ch], "_", i, "_", marker.allele.text, sep = "")
						num.blocks <- num.blocks + 1
						}
					}
				}
			}
		}			
	
	if(collapse.linked.markers){
		data.obj$linkage.blocks.collapsed <- linkage.blocks
		}else{
		data.obj$linkage.blocks.full <- linkage.blocks
		}
	return(data.obj)
	
}