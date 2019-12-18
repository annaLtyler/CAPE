#this function plots chromosome bars with labels for plotNetworkDO
# start = 628; start = 5
plot.chr <- function(chr.rad, label.rad, chr, start, pts.per.chr, gap){
	
	
		chr.circ <- get.circle(chr.rad)
		label.circ <- get.circle(label.rad)

		# plot.new()
		# plot.window(xlim = c(-4, 4), ylim = c(-4, 4))	
		# #test the chr circ and label circ
		# points(chr.circ$x, chr.circ$y)
		# points(label.circ$x, label.circ$y, col = "gray")
	
		for(i in 1:length(chr)){
			chr.x.coord <- chr.circ$x[start:(start+pts.per.chr[i]-1)]
			chr.y.coord <- chr.circ$y[start:(start+pts.per.chr[i]-1)]
			
			points(chr.x.coord, chr.y.coord, type = "l", lwd = main.lwd)

			label.x.coord <- label.circ$x[start:(start+pts.per.chr[i]-1)]
			label.y.coord <- label.circ$y[start:(start+pts.per.chr[i]-1)]
			
			label.coord <- get.arc.mid(label.x.coord, label.y.coord)
			
			if(names(chr)[i] == "covar"){
				text.adj = 0
				}else{
				text.adj = 1
				}
				
			text(label.coord[1], label.coord[2], chr.names[i], adj = text.adj, cex = label.cex)
			
			start = start + pts.per.chr[i] + gap
			}
	
	
	
}