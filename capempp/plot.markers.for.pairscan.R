#This function plots the lod score from a singlescan
#with the markers selected for the pairscan marked

plot.markers.for.pairscan <- function(data.obj, singlescan.obj, by.LOD = TRUE, min.lod = 0, smoothness = 0.02, color.scheme = c("DO/CC", "other")){

	plot.peaks = TRUE
	# library(RColorBrewer)
	# cols <- brewer.pal(12, "Set3")
	cols <- c("grey", "white")

	markers <- rownames(singlescan.obj$locus.lod.scores)
	chr <- get.marker.chr(data.obj, markers)
	chr <- chr[which(chr != 0)]
	u_chr <- unique(chr)
	num.pheno <- ncol(singlescan.obj$locus.lod.scores)
	all.bins <- vector(mode = "list", length = num.pheno)
	names(all.bins) <- colnames(singlescan.obj$locus.lod.scores)

	get.marker.max <- function(allele.row){
		return(apply(allele.row, 1, function(x) max(abs(x))))
		}

	get.max.allele <- function(allele.row){
		return(colnames(allele.row)[apply(allele.row, 1, function(x) which.max(abs(x)))])
		}

	if(by.LOD){
		max.y <- max(singlescan.obj$singlescan.t.stats)
		curve.by.pheno <- t(apply(singlescan.obj$singlescan.t.stats, 1, get.marker.max))
		}else{
		max.y <- max(singlescan.obj$locus.lod.scores)
		curve.by.pheno <- singlescan.obj$locus.lod.scores
		}
	max.allele <- t(apply(singlescan.obj$singlescan.t.stats, 1, get.max.allele))
	
	
	for(l in 1:num.pheno){
		plot.new()
		plot.window(xlim = c(0, length(markers)), ylim = c(0, max.y))

		chr.bins <- vector(mode = "list", length = length(u_chr))
		names(chr.bins) <- u_chr

		for(ch in 1:length(u_chr)){
			
			chr.locale <- which(chr == u_chr[ch])
			lod <- curve.by.pheno[chr.locale,l]
			chr.alleles <- max.allele[chr.locale,l]
			
			smoothed <- loess.smooth(c(1:length(lod)), lod, span = smoothness)
			if(plot.peaks){
				points(chr.locale, lod, type = "l")
				}
			
			smoothed.y <- c(0, smoothed$y, 0)	
			smoothed.left <- smoothed$y - head(smoothed.y, length(smoothed$y))
			smoothed.right <- smoothed$y - tail(smoothed.y, length(smoothed$y))
		
			peak.locale <- intersect(which(smoothed.left > 0), which(smoothed.right > 0))
			peak.x <- smoothed$x[peak.locale]
			peak.y <- smoothed$y[peak.locale]	
			
			trough.locale <- intersect(which(smoothed.left < 0), which(smoothed.right < 0))
			trough.x <- smoothed$x[trough.locale]
			trough.y <- smoothed$y[trough.locale]
			
			padded.trough <- unique(c(1, trough.locale, length(smoothed$y)))
			bins <- vector(mode = "list", length = (length(padded.trough)-1))
			for(j in 1:(length(padded.trough)-1)){
				poly.x <- round(smoothed$x[padded.trough[j]]:smoothed$x[padded.trough[j+1]])
				poly.y <- lod[poly.x]
				poly.y[1] <- 0;poly.y[length(poly.y)] <- 0
				bins[[j]] <- cbind(markers[chr.locale[poly.x]], lod[poly.x], chr.alleles[poly.x])
				if(plot.peaks){
				polygon(x = chr.locale[poly.x], y = poly.y, col = cols[j%%length(cols)])
				}
				}#end looping through bins
				chr.bins[[ch]] <- bins
				par(xpd = TRUE)
				text(x = mean(chr.locale), y = max.y*-0.05, labels = u_chr[ch])
				par(xpd = FALSE)
				abline(v = max(chr.locale), lwd = 3, col = "gray", lty = 2)
			}#end looping through chromosomes
			axis(2)
			mtext(colnames(curve.by.pheno)[l])
			mtext("LOD score", side = 2, line = 2.5)
			all.bins[[l]] <- chr.bins
	
		selected.markers <- colnames(data.obj$geno.for.pairscan)
		#get the x coordinates for these markers
		split.markers <- strsplit(selected.markers, "_")
		just.markers <- unlist(lapply(split.markers, function(x) x[1]))
		just.alleles <- unlist(lapply(split.markers, function(x) x[2]))
	
		alleles <- unique(just.alleles)
		x.locale <- match(just.markers, markers)
		allele.cols <- get.allele.colors(color.scheme, alleles)
		allele.idx <- unlist(lapply(just.alleles, function(x) which(alleles == x)))
		all.cols <- allele.cols[allele.idx,3]
		points(x = x.locale, y = rep(0, length(x.locale)), pch = "|", col = all.cols)
		abline(h = min.lod)
		}#end looping through traits
	
	
	
	
}