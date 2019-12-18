	#this function returns the points on a line between two points
	get.line <- function(x0, y0, x1, y1, dens = circle.dens){
		x.pts <- segment.region(x0, x1, 1/dens, alignment = "ends")
		slope <- (y1-y0)/(x1-x0)
		y.pts <- slope*(x.pts - x0) + y0
		result <- list(x.pts, y.pts); names(result) <- c("x", "y")
		return(result)
		}
		
