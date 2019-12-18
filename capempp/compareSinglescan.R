#This script makes effects plots for multi-allelic 
#1D scans.
#The view can either be condensed to show an overview
#of the full genome, or spread out to show a detailed
#view of each chromosome. The default view is "overview"

compareSinglescan <- function(data.obj, geno.obj = NULL, singlescan.obj, chr = NULL, traits = NULL, standardized = TRUE, allele.labels = NULL, include.covars = TRUE, view = c("overview", "detailed"), above.t = NULL){
	
	if(is.null(chr)){
		chr <- sort(unique(data.obj$chromosome))
		if(!include.covars && length(which(chr == "0")) > 0){
			chr <- chr[-(which(chr == "0"))]
			}
		}
	
	if(is.null(above.t)){
		above.t <- 0
		}

	geno <- get.geno(data.obj, geno.obj)

	#Get the dimension names to minimize confusion	
	mouse.dim <- which(names(dimnames(geno)) == "mouse")
	locus.dim <- which(names(dimnames(geno)) == "locus")
	allele.dim <- which(names(dimnames(geno)) == "allele")


	all.chromosomes <- data.obj$chromosome
		
	
	if(!standardized){
		results <- singlescan.obj$singlescan.effects
		plot.type.label <- "beta"
		}else{
			results <- singlescan.obj$singlescan.t.stats
			plot.type.label <- "t.stat"
			}
	
	
	if(is.null(traits)){
		traits <- dimnames(results)[[2]]
		}
	results.el <- which(dimnames(results)[[2]] %in% traits)
	
	if(length(results.el) < length(traits)){
		if(length(results.el) > 0){
			not.found <- traits[-results.el]
			}else{
			not.found <- traits
			}
		message("I couldn't find the following traits:")
		cat(not.found, sep = "\n")
		return()
		}
	
	#subset the results based on which chromosomes
	#are desired.
	chr.locale <- which(all.chromosomes %in% chr)
	sub.results <- results[chr.locale,,,drop=FALSE]

	#subset the results based on the effect size threshold
	high.locale <- which(abs(sub.results) >= above.t, arr.ind = TRUE)
	# print(str(high.locale))
	u_high <- sort(unique(high.locale[,1]))
	
	high.results <- array(NA, dim = c(length(u_high), dim(sub.results)[[2]], dim(sub.results)[[3]]))
	if(length(high.results) == 0){stop("There are no markers above this threshold.")}
	
	all.rows <- rep(NA, length(u_high))
	for(i in 1:length(u_high)){
		m.locale <- which(high.locale[,1] == u_high[i])
		for(j in 1:length(m.locale)){
			high.results[i,high.locale[m.locale[j],2], high.locale[m.locale[j],3]] <- sub.results[high.locale[m.locale[j],1], high.locale[m.locale[j],2], high.locale[m.locale[j],3]]
			all.rows[i] <- rownames(sub.results)[high.locale[m.locale[j]]]
			}
		}
	dimnames(high.results)[[1]] <- all.rows
	dimnames(high.results)[[2]] <- dimnames(sub.results)[[2]]
	dimnames(high.results)[[3]] <- dimnames(sub.results)[[3]]
	
	sub.results <- high.results
	
	if(length(results) == 0){
		stop("You must run singlescan.R before plotting effects.")
		}

	if(length(dim(results)) < 3){
		stop("This function is for plotting effects of multiple alleles.\nYou only have two alleles at each locus.")
		}
	
	#For each phenotype, we want to plot the effects of the
	#presence of each parent allele across the genome in its
	#own color
	
	num.loci <- dim(sub.results)[[1]]
	num.alleles <- dim(sub.results)[[3]]
	alleles <- dimnames(sub.results)[[3]]
	allele.colors <- get.allele.colors(color.scheme, alleles)
	used.alleles <- allele.colors[,2]
	allele.labels <- allele.colors[,1]
	
	
	phenos.scanned  <- dimnames(results)[[2]]
	
	
	
	if(length(grep("over", view)) > 0){
		view <- "overview"
		mark.what = NULL
		}else{
		view <- "detailed"
		mark.what = NULL
		}
	
		
	if(view == "overview"){
		layout.mat <- matrix(c(1:(length(results.el)*2)), ncol = 1)
		layout(layout.mat, heights = rep(c(0.85, 0.15), length(results.el)))
		# layout.show(length(results.el)*2)
		
		for(i in results.el){
			# dev.new(width = 15, height = 5)
			if(plot.type.label == "t.stat"){
				par(mar = c(0,5,5,2))
				plot.new()
				# plot.window(xlim = c(1,num.loci), ylim = c(min(abs(sub.results), na.rm = TRUE), max(abs(sub.results), na.rm = TRUE)))
				plot.window(xlim = c(1,num.loci), ylim = c(0, max(abs(sub.results), na.rm = TRUE)))
				}else{
					par(mar = c(0,5,5,2))
					plot.new()
					plot.window(xlim = c(1,num.loci), ylim = c(min(sub.results, na.rm = TRUE), max(sub.results, na.rm = TRUE)))
					}
				
			for(j in 1:num.alleles){
				#pull out the effects of the presence of
				#allele j on phenotype i
				allele.effects <- as.vector(sub.results[,i,j])
				if(plot.type.label == "t.stat"){ #plot the absolute value of the t.statistics
					points(abs(allele.effects), col = allele.colors[j,3], type = "l", lwd = 1)
					}else{
						points(allele.effects, col = allele.colors[j,3], type = "l", lwd = 1)	
						}
						
				if(plot.type.label == "t.stat"){
					lines(x = c(1,num.loci), y = rep(data.obj$pairscan.thresh, 2), lty = 1, col = "darkgray")
					lines(x = c(1,num.loci), y = rep(data.obj$covar.thresh, 2), lty = 2, col = "darkgray")
					# abline(h = data.obj$covar.thresh, lty = 2, col = "darkgray")
					par(xpd = TRUE)
					text(x = num.loci*1.01, y = data.obj$pairscan.thresh, labels = paste("p =", data.obj$alpha.for.pairs), cex = 0.5, adj = 0)
					text(x = num.loci*1.01, y = data.obj$covar.thresh, labels = paste("p =", data.obj$alpha.for.covar), cex = 0.5, adj = 0)
					par(xpd = FALSE)
					}
				} #end looping over alleles
			
			abline(h = 0)
			axis(2)
			# axis(1, labels = FALSE)
			mtext(paste("Effect Relative to Allele", ref.allele), side = 2, line = 2.5)
			mtext(phenos.scanned[i], line = 2.7)
			
			if(plot.type.label == "t.stat"){
			legend(0, (max(sub.results, na.rm = TRUE)*1.2), legend = used.alleles, col = allele.colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
			}
			
			#put in lines for chromosome boundaries
			for(ch in 1:length(chr)){
				abline(v = max(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr)), lty = 3)
				}	
		
			if(!is.null(mark.what)){
				points(ind.locale, rep(max(sub.results)*1.02, length(ind.locale)), col = "red", pch = "*")
				}
	
			par(mar = c(0,5,0,2))	
			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(0,1))		
			#indicate where chromosomes are and rewrite the
			#phenotype for each, so we can see it on really
			#zoomed in plots
			for(ch in 1:length(chr)){
				#find the mean position of where each chromosome starts and stops
				mid.pos <- mean(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr))
				text(mid.pos, 0.7, labels = chr[ch], xpd = TRUE, cex = 1.5)
				}
				
			mtext("Chromosome", side = 1, line = -1.6, cex = 1.2)	
							
			} #end looping over phenotypes
		
		
		}else{ #switch cases from overview to detailed view
	
	

		layout.mat <- matrix(1:(length(results.el)*3), ncol = 1)
		layout(layout.mat, heights = rep(c(0.07, 0.78, 0.15), length(results.el)))
		# layout.show(length(results.el)*3)		
		
		for(i in results.el){
			# dev.new(width = length(chr)*10, height = 5)
			#plot the phenotype name above each chromosome, so
			#we can see it on really zoomed in plots
			par(mar = c(0,5,0,2))	
			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(0,1))		
			#indicate where chromosomes are and rewrite the
			#phenotype for each, so we can see it on really
			#zoomed in plots
			for(ch in 1:length(chr)){
				#find the mean position of where each chromosome starts and stops
				mid.pos <- mean(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr))
				text(mid.pos, 0.5, labels = phenos.scanned[i], xpd = TRUE, cex = 2)
				}	


			par(mar = c(0,5,5,2))
			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(min(sub.results, na.rm = TRUE), max(sub.results, na.rm = TRUE)))
				
				for(j in 1:num.alleles){
					#pull out the effects of the presence of
					#allele j on phenotype i
					allele.effects <- as.vector(sub.results[,i,j])
					if(plot.type.label == "t.stat"){ #plot the absolute value of the t.statistics
						points(abs(allele.effects), col = allele.colors[j,3], type = "l", lwd = 2)
						}else{
							points(allele.effects, col = allele.colors[j,3], type = "l", lwd = 2)	
							}
					}

			abline(h = 0)
			axis(2)
			# axis(1, labels = FALSE)
			mtext(paste("Effect Relative to Allele", ref.allele), side = 2, line = 2.5, cex = 1.5)
			# mtext(phenos.scanned[i], line = 1.5, cex = 2)
			legend(0, (max(sub.results, na.rm = TRUE)*1.2), legend = used.alleles, col = allele.colors[,3], lty = 1, lwd = 3, xpd = TRUE, horiz = TRUE)
			#put in lines for chromosome boundaries
			for(ch in 1:length(chr)){
				abline(v = max(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr)), lty = 3)
				}	
	
			if(plot.type.label == "t.stat"){
				abline(h = data.obj$pairscan.thresh, lty = 1, col = "darkgray")
				abline(h = data.obj$covar.thresh, lty = 2, col = "darkgray")

				}
	
		
			if(!is.null(mark.what)){
				points(ind.locale, rep(max(sub.results)*1.02, length(ind.locale)), col = "red", pch = "*")
				}

	
			par(mar = c(0,5,0,2))	
			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(0,1))		
			#indicate where chromosomes are and rewrite the
			#phenotype for each, so we can see it on really
			#zoomed in plots
			for(ch in 1:length(chr)){
				#find the mean position of where each chromosome starts and stops
				mid.pos <- mean(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr))
				text(mid.pos, 0.5, labels = paste("Chr", chr[ch]), xpd = TRUE, cex = 2)
				}	
		
			}

		}
		
	
	}