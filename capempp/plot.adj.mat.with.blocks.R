
plot.adj.mat.with.blocks <- function(adj.mat, l.blocks, plot.label = NULL, color.breaks = c(0.5)){						
	#l.blocks from diagonal means assign the diagonals
	#not the columns/rows to blocks. First we need to
	#calculate the dimensions of the block from the 
	#diagnoal coordinates
	#there are (n*2)-1 diagonal lines in the matrix
	
	# my.palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")

	# image(1:dim(adj.mat)[1], 1:dim(adj.mat)[2], adj.mat, xlim = c(0,(dim(adj.mat)[1]+1)), ylim = c(0,(dim(adj.mat)[1]+1)), col = my.palette(50), axes = FALSE, ylab = "", xlab = "", main = plot.label)
	color.order <- c("blue", "green", "yellow", "orange", "gray")
	
	layout(matrix(c(1,2), nrow = 1), widths = c(1, 0.2))
	if(is.null(color.breaks)){
	imageWithText(mat = rotate.mat(rotate.mat(rotate.mat(adj.mat))), show.text = FALSE, col.scale = color.order, global.color.scale = TRUE)
		}else{
	imageWithText(rotate.mat(rotate.mat(rotate.mat(adj.mat))), show.text = FALSE, split.at.vals = TRUE, split.points = color.breaks, col.scale = color.order, global.color.scale = TRUE)
		}
	
	start.block = 0.5
	block.locale <- sapply(1:max(l.blocks), function(x) which(l.blocks == x))

	#outline each block
	for(b in 1:length(block.locale)){
		
		end.block <- start.block + (length(block.locale[[b]])/2)
		segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
		segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
		segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
		segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
		start.block <- end.block
		} #end outlining blocks


	if(is.null(color.breaks)){
		imageWithTextColorbar(rotate.mat(rotate.mat(rotate.mat(adj.mat))), col.scale = color.order, cex = 1.2, global.color.scale = TRUE)
		}else{
		imageWithTextColorbar(rotate.mat(rotate.mat(rotate.mat(adj.mat))), split.at.vals = TRUE, split.points = color.breaks, col.scale = color.order, cex = 1.2, global.color.scale = TRUE)
		}


	}
