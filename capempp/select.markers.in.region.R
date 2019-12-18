#This function selects markers in a region
#based on the singlescan data. This is 
#for going from a significant marker to
#a region for gene hunting as an alternative
#to collapsing linkage blocks
#anchor.marker = "2:165241741_D"; allele = "D";trait = "me3"

select.markers.in.region <- function(data.obj, singlescan.obj, anchor.marker, trait, allele, window.size = rep(100, 3), plot.peaks = FALSE, verbose = TRUE){


	cols <- c("grey", "white")

	if(class(singlescan.obj) == "list"){ 
		results <- singlescan.obj$singlescan.t.stats #an actual singlescan object
		}else{
		results <- singlescan.obj #a singlscan matrix for calculating pairscan null distribution
		}


	
	covar.info <- get.covar(data.obj)
	results.no.covar <- abs(results[which(!rownames(results) %in% covar.info$covar.names),,])
	
	
	marker.name = strsplit(anchor.marker, "_")[[1]][1]
	marker.locale <- which(dimnames(results.no.covar)[[1]] == marker.name) 
	marker.chr <- get.marker.chr(data.obj, marker.name)
	all.chr <- get.marker.chr(data.obj, dimnames(results.no.covar)[[1]])
	just.marker.chr.locale <- which(all.chr == marker.chr)
	
	trait.locale <- which(dimnames(results.no.covar)[[2]] == trait)
	strain.locale <- which(dimnames(results.no.covar)[[3]] == allele)
	just.trait.curve <- results.no.covar[just.marker.chr.locale,trait.locale,strain.locale]
	# plot(just.trait.curve, type = "l")
	#===============================================================
	#internal functions
	#===============================================================
		
		#how may peaks are above a given cutoff?
		num.peaks <- function(allele.curves, bins, cutoff){
			filtered.bins <- bins
			filtered.bins[which(allele.curves < cutoff)] <- NA
			num.peaks <- apply(filtered.bins, 2, function(x) length(unique(x))-1)
			return(num.peaks)
			}
			
		#how many alleles are above a given t stat cutoff
		num.markers <- function(allele.curves, cutoff){
			filtered.curves <- allele.curves
			filtered.curves[which(allele.curves < cutoff)] <- NA
			num.markers <- apply(filtered.curves, 2, function(x) length(which(!is.na(x))))
			return(num.markers)
			}

		#how many markers are in each peak at a given cutoff?
		markers.per.peak <- function(allele.curves, bins, cutoff){
			result.mat <- matrix(0, ncol = ncol(allele.curves), nrow = max(bins, na.rm = TRUE))
			rownames(result.mat) <- 1:nrow(result.mat)
			colnames(result.mat) <- colnames(allele.curves)
			filtered.bins <- bins
			filtered.bins[which(allele.curves < cutoff)] <- NA
			for(i in 1:ncol(filtered.bins)){
				counts <- table(filtered.bins[,i])
				result.mat[names(counts),i] <- counts
				}
			return(result.mat)
			}
	
		sample.peaks <- function(pheno.results, num.per.peak, bins){
			sampled.markers <- vector(mode = "list", length = ncol(pheno.results))
			names(sampled.markers) <- colnames(pheno.results)
			for(i in 1:ncol(pheno.results)){
				#figure out which peaks we will sample from
				allele.markers <- NULL
				peaks.which <- which(num.per.peak[,i] > 0)
				if(length(peaks.which) > 0){
					for(j in 1:length(peaks.which)){
						#in each peak, pick the max, and sample the rest
						marker.locale <- which(bins[,i] == peaks.which[j])
						allele.markers <- c(allele.markers, marker.locale[which.max(pheno.results[marker.locale,i])])
						num.to.sample <- num.per.peak[peaks.which[j],i] - 1 #take off the maximum marker
						if(num.to.sample > 0){ #if there are still markers to get after grabbing the max
							unif.markers <- round(runif(num.to.sample, min = min(marker.locale), max = max(marker.locale)))
							allele.markers <- c(allele.markers, unif.markers)
							}#end case for sampling peak uniformly
						}#end looping through peaks for one allele
					sampled.markers[[i]] <- rownames(pheno.results)[sort(allele.markers)]
					}#end looping through alleles
				}
			return(sampled.markers)
			}
	
		#===============================================================

	bins <- bin.curve(just.trait.curve, plot.peaks = plot.peaks, window.size = window.size)$bins

	new.marker.locale <- which(dimnames(results.no.covar)[[1]][just.marker.chr.locale] == marker.name)
	marker.bin <- bins[new.marker.locale]
	other.marker.locale <- which(bins == marker.bin)
	other.markers.in.bin <- dimnames(results.no.covar)[[1]][just.marker.chr.locale[other.marker.locale]]
	if(plot.peaks){
		points(x = other.marker.locale, y = rep(0, length(other.marker.locale)), pch = "*", col = "red")
		abline(v = new.marker.locale)
		}
	return(other.markers.in.bin)
	
}