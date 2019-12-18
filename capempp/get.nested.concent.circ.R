#This function makes a named list of nested 
#concentric circles such that each major 
#circle consists of a set of minor circles
#this is used for plotting allele contributions
#to traits in plotNetworkDO

	get.nested.concent.circ <- function(major.list, minor.list, start.rad = 1.05, major.gap = 0.05, minor.gap = 0.01){

		trait.circ <- vector(mode = "list", length = length(major.list))
		names(trait.circ) <- major.list
		num.minor = length(minor.list)
		
		for(i in 1:length(major.list)){
			for(j in 1:num.minor){
				effect.circ = get.circle(start.rad)
				trait.circ[[i]][[j]] <- effect.circ
				start.rad <- start.rad + minor.gap
				}
			names(trait.circ)[[i]] <- minor.list[i]
			start.rad <- start.rad + major.gap
			}
			

	# #test plot
	# allele.colors <- c("gold", "black", "salmon", "blue", "lightblue", "green", "red", "purple")
	# plot.new()
	# plot.window(xlim = c((1+start.rad), (1-start.rad)), ylim = c((1+start.rad), (1-start.rad)))
	# for(i in 1:length(trait.circ)){
		# for(j in 1:length(minor.list)){
			# ind.col <- colorRampPalette(c("white", allele.colors[j]))
			# rnd.vals <- rpois(length(trait.circ[[i]][[j]]$x), lambda = 1)/10
			# colorRamp <- ind.col(length(rnd.vals))
			# points(trait.circ[[i]][[j]]$x, trait.circ[[i]][[j]]$y, col = colorRamp, cex = 0.1)
			# }
		# }

		return(trait.circ)


		}


