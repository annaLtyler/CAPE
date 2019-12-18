#This function plots the trait circles which have gaps at
#the right-hand sice with labels

	plot.nested.trait.circ <- function(trait.circ, label.gap){
		
		all.min <- min(unlist(trait.circ))
		all.max <- max(unlist(trait.circ))
		
		num.pheno = length(trait.circ)
		
		pheno.label.starts <- round(segment.region(1, label.gap, num.points = num.pheno, alignment = "center"))
		
		#put the labels between the outermost trait circle and the edge of the plot
		num.traits <- length(trait.circ)
		num.nested <- lapply(trait.circ, length)
		max.end <- as.numeric(tail(num.nested, 1))
		label.x <- max(trait.circ[[num.traits]][[max.end]]$x) + (plot.dim[2] - max(trait.circ[[num.traits]][[max.end]]$x))*0.25
		arrow.start.x <- max(trait.circ[[num.traits]][[max.end]]$x) + (plot.dim[2] - max(trait.circ[[num.traits]][[max.end]]$x))*0.2
	
		# plot.new()
		# plot.window(xlim = c(all.min, all.max*1.5), ylim = c(all.min, all.max))
		for(tr in length(trait.circ):1){
			#add light gray bars to help the eye track the phenotypes
			#they should be staggered so the label sticks don't overlap
			#the allele colors will be plotted over the top of the gray bars
			for(n in 1:length(trait.circ[[tr]])){
				circ.x <- trait.circ[[tr]][[n]]$x[pheno.label.starts[tr]:length(trait.circ[[tr]][[n]]$x)]
				circ.y <- trait.circ[[tr]][[n]]$y[pheno.label.starts[tr]:length(trait.circ[[tr]][[n]]$y)]
				points(circ.x, circ.y, type = "l", col = "lightgray", lwd = main.lwd)
				}
			
			#and add phenotype labels
			arrow.end.x <- trait.circ[[tr]][[n]]$x[pheno.label.starts[tr]]
			arrow.y <- trait.circ[[tr]][[n]]$y[pheno.label.starts[tr]]
			segments(x0 = arrow.start.x, x1 = arrow.end.x, y0 = arrow.y, lwd = main.lwd, col = "lightgray")
			text(label.x, arrow.y, labels = names(trait.circ)[tr], adj = 0)
			}
		
	}
