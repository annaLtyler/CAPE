#This function plots a subnetwork 
#related to a given marker
# marker.name = "Sex"; p.or.q = 0.05;  collapsed.net = TRUE; pos.col = "brown"; neg.col = "blue"; bg.col = "gray"; light.dark = "l"; node.border.lwd = 1; layout.matrix = NULL; zoom = 1; xshift = 0; yshift = 0; node.radius = 1; label.nodes = TRUE; label.offset = 0; label.cex = 1; legend.radius = 1; legend.cex = 1; legend.position = "topleft"; arrow.offset = node.radius; arrow.length = 0.2; edge.lwd = 2

plotNetwork.subnet <- function(data.obj, marker.name, p.or.q = 0.05,  collapsed.net = TRUE, pos.col = "brown", neg.col = "blue", bg.col = "gray", light.dark = "l", node.border.lwd = 1, layout.matrix = NULL, zoom = 1, xshift = 0, yshift = 0, node.radius = 1, label.nodes = TRUE, label.offset = 0, label.cex = 1, legend.radius = 1, legend.cex = 1, legend.position = "topleft", arrow.offset = node.radius, arrow.length = 0.2, edge.lwd = 2){

	new.data <- data.obj

	if(collapsed.net){
		net <- data.obj$collapsed.net
		}else{
		net <- data.obj$full.net	
		}

	just.int <- net[,1:nrow(net)]
	just.main <- net[,(nrow(net)+1):ncol(net)]

	marker.locale <- which(rownames(net) == marker.name)


	#plot networks by phenotype
	if(length(marker.locale) > 0){
		int.locale <- sort(unique(c(which(just.int[marker.locale,] != 0), which(just.int[,marker.locale] != 0), marker.locale)))

			new.net <- cbind(just.int[int.locale,int.locale,drop=FALSE], just.main[int.locale,,drop=FALSE])
			if(collapsed.net){
				new.data$collapsed.net <- new.net
				}else{
				new.data$full.net <- new.net	
				}
			plotNetwork2(new.data, main = colnames(cross$pheno)[i], p.or.q = p.or.q,  collapsed.net = collapsed.net, pos.col = pos.col, neg.col = neg.col, bg.col = bg.col, light.dark = light.dark, node.border.lwd = node.border.lwd, layout.matrix = layout.matrix, zoom = zoom, xshift = xshift, yshift = yshift, node.radius = node.radius, label.nodes = label.nodes, label.offset = label.offset, label.cex = label.cex, legend.radius = legend.radius, legend.cex = legend.cex, legend.position = legend.position, arrow.offset = arrow.offset, arrow.length = arrow.length, edge.lwd = edge.lwd)
			}
		invisible(new.net)

}