#This function gets the middle point of an arc
#just by getting the middle of the x and y 
#coords. it is used in plotNetwork

get.arc.mid <- function(arc.x, arc.y){
	
	mid.point <- round(length(arc.x)/2)
	result <- c("x" = arc.x[mid.point], "y" = arc.y[mid.point])
	return(result)
	
	
}