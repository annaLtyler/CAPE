#This script makes effects plots for multi-allelic 
#1D scans.
#The view can either be condensed to show an overview
#of the full genome, or spread out to show a detailed
#view of each chromosome. The default view is "overview"

plotSinglescanRaster <- function(data.obj, geno.obj = NULL, singlescan.obj, chr = NULL, traits = NULL, standardized = TRUE, allele.labels = NULL, include.covars = TRUE, color.thresh = NULL){
	
	if(is.null(chr)){
		chr <- sort(unique(data.obj$chromosome))
		if(!include.covars && length(which(chr == "0")) > 0){
			chr <- chr[-(which(chr == "0"))]
			}
		}
	

	geno <- get.geno(data.obj, geno.obj)

	#Get the dimension names to minimize confusion	
	mouse.dim <- which(names(dimnames(geno)) == "mouse")
	locus.dim <- which(names(dimnames(geno)) == "locus")
	allele.dim <- which(names(dimnames(geno)) == "allele")


	all.chromosomes <- data.obj$chromosome
	lod.scores <- singlescan.obj$locus.lod.scores
	
		
	
	if(!standardized){
		results <- singlescan.obj$singlescan.effects
		plot.type.label <- "beta"
		}else{
			results <- abs(singlescan.obj$singlescan.t.stats)
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
	sub.results <- results[chr.locale,,]
	lod.scores <- lod.scores[chr.locale,]


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
	ref.allele <- singlescan.obj$ref.allele
	alleles <- dimnames(geno)[[allele.dim]]
	ref.allele.locale <- which(alleles == ref.allele)
	# allele.colors <- c("gold", "black", "salmon", "blue", "lightblue", "green", "red", "purple")
	# used.colors <- allele.colors[-ref.allele.locale] #take out the color of the ref.allele
	
	#the color ramp start and stop for each of the colors
	#higher colors are lighter
	color.ramps <- list(c(rgb(255,215,0, max = 255), rgb(255,253,222, max = 255)), #gold
						c(rgb(1,1,1, max = 255), rgb(243,243,243, max = 255)), #black
						c(rgb(255,128,51, max = 255), rgb(255,204,229, max = 255)), #orange instead of salmon
						c(rgb(0,0,153, max = 255), rgb(204,204,255, max = 255)), #dark blue
						c(rgb(51,153,255, max = 255), rgb(204,229,255, max = 255)), #light blue
						c(rgb(0,204,0, max = 255), rgb(204,255,229, max = 255)), #green
						c(rgb(255,0,0, max = 255), rgb(255,204,204, max = 255)), #red
						c(rgb(127,0,255, max = 255), rgb(229,204,255, max = 255))) #purple

	#lower colors are lighter
	# color.ramps <- list(c(rgb(255,253,222, max = 255), rgb(255,247,0, max = 255)), #gold
						# c(rgb(243,243,243, max = 255), rgb(1,1,1, max = 255)), #black
						# c(rgb(255,204,229, max = 255), rgb(255,128,51, max = 255)), #orange instead of salmon
						# c(rgb(204,204,255, max = 255), rgb(0,0,153, max = 255)), #dark blue
						# c(rgb(204,229,255, max = 255), rgb(51,153,255, max = 255)), #light blue
						# c(rgb(204,255,229, max = 255), rgb(0,204,0, max = 255)), #green
						# c(rgb(255,204,204, max = 255), rgb(255,0,0, max = 255)), #red
						# c(rgb(229,204,255, max = 255), rgb(127,0,255, max = 255))) #purple

	
	color.ramps[[ref.allele.locale]] <- NULL
	if(is.null(allele.labels)){
		used.alleles <- alleles[-ref.allele.locale]
		}else{
			if(length(allele.labels) == dim(geno)[allele.dim]){
				used.alleles <- allele.labels[-ref.allele.locale]
				}else{
				used.alleles <- allele.labels	
				}
			}
	
	
	phenos.scanned  <- dimnames(results)[[2]]
	
	view <- "overview"
		
	#plot each allele in a different window, so we don't have to 
	#worry about getting everything smooshed together
	#we need one window for the title + one for each allele + a 
	#window for chromosome labels and a window for the word "chromosome"
	layout.mat <- matrix(1:(length(color.ramps)+3), ncol = 1)
	

		for(i in results.el){
			# dev.new(width = 15, height = 5)
	
			# quartz(width = 10, height = 3)
			layout(layout.mat)
			# layout.show(dim(layout.mat)[1])

			#margins for the windows with the alleles
			par(mar = c(0,5,0,2))			

			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(0,1))
			text(x = num.loci/2, y = 0.5, labels = phenos.scanned[i], cex = 2)

			par(xpd = TRUE)
			for(j in 1:num.alleles){
				plot.new()
				plot.window(xlim = c(1,num.loci), ylim = c(0.45,0.55))
				plot.dim <- par("usr")
				#pull out the effects of the presence of
				#allele j on phenotype i
				#and build a matrix of colors
				allele.effects <- as.vector(sub.results[,i,j])
				if(plot.type.label == "beta"){
					allele.effects <- as.vector(sub.results[,i,j]) + min(as.vector(sub.results[,i,j])) #move everything above 0
					}
				above.thresh.locale <- which(allele.effects >= color.thresh)
				allele.effects[above.thresh.locale] <- min(allele.effects) #first set these to a low value
				col.ramp <- colorRamp(colors = color.ramps[[j]])
				allele.rgb <- col.ramp((allele.effects/max(allele.effects))^(1/2))
				# allele.rgb <- col.ramp(allele.effects/max(allele.effects))
				allele.col <- apply(allele.rgb, 1, function(x) rgb(x[1], x[2], x[3], max = 255))
				
				#change the color of the values above the threshold to white
				allele.col[above.thresh.locale] <- rgb(255, 255, 255, max = 255)
				
				# points(x = 1:length(allele.effects), y = rep(num.alleles-j+1, length(allele.effects)), col = allele.col, type = "p", pch = "|", cex = 1)
				# text(x = plot.dim[1], y = num.alleles-j+1, labels = used.alleles[j])	
				# points(x = 1:length(allele.effects), y = rep(0.5, length(allele.effects)), col = allele.col, type = "p", pch = "|", cex = 1)
				segments(x0 = 1:length(allele.effects), y0 = rep(0.45, length(allele.effects)), y1 = rep(0.55, length(allele.effects)), col = allele.col)
				text(x = plot.dim[1], y = 0.5, labels = used.alleles[j])	

				#put in lines for chromosome boundaries
				abline(v = 0)
				for(ch in 1:length(chr)){
					line.type = 1
					abline(v = max(which(all.chromosomes == chr[ch])) - min(which(all.chromosomes %in% chr))+1, lty = line.type)
					}	

				#also put lines at the top and bottom
				if(j == 1){
					segments(x0 = 1, y0 = 0.55, x1 = length(allele.effects), col = "black", lwd = 1, lty = 1)
					}
				if(j == num.alleles){
					segments(x0 = 1, y0 = 0.45, x1 = length(allele.effects), col = "black", lwd = 1, lty = 1)
					}


				} #end looping over alleles
			# par(xpd = FALSE)
		
		
			#margins for the window with the chromosome labels
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
			
			plot.new()
			plot.window(xlim = c(1,num.loci), ylim = c(0,1))
			text(x = num.loci/2, y = 0.5, labels = "Chromosome", cex = 2)	
			} #end looping over phenotypes
		
	
	}