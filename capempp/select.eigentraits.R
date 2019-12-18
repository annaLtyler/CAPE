#This script is used to select individual eigentraits
#after viewing the results of the svd
#The script defaults to eigentraits 1 and 2

select.eigentraits <- function(data.obj, traits.which = c(1,2)){
	
	check.underscore(data.obj)
	# check.bad.markers(data.obj)
	
	ET <- data.obj$ET
	selected.ET <- ET[,traits.which]
	data.obj$ET <- selected.ET
	return(data.obj)
	
}