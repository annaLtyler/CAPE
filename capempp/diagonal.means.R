#This function finds the mean of the diagonals of a matrix
#it sweeps through a matrix finding the means of each sucessive
#diagonal.
#the function starts either in the top left or top right of the 
#matrix


diagonal.means <- function(mat, start.pt = c("tl", "br", "tr", "bl")){
	
	if(nrow(mat) != ncol(mat)){
		stop('This is not a square matrix')
		}
		
	if(length(start.pt) > 1){
		start.pt <- "tl"
		}
				
	all.means <- rep(NA, (ncol(mat)*2)-1)
	num.row <- nrow(mat)
	mat1 <- matrix(1:num.row, num.row, num.row, byrow = TRUE)
	mat2 <- matrix(1:num.row, num.row, num.row, byrow = FALSE)

	if(start.pt == "bl" || start.pt == "tr"){

		
		diff.mat <- mat1-mat2
		min.idx <- min(diff.mat)
		max.idx <- max(diff.mat)
		all.idx <- min.idx:max.idx
		
		for(i in 1:length(all.idx)){
			report.progress(i, length(all.idx))
			idx.locale <- which(diff.mat == all.idx[i])
			diag.v <- mat[idx.locale]
			all.means[i] <- mean(diag.v, na.rm = TRUE)
			}

		if(start.pt == "tl"){
			return(rev(all.means))
			}
		if(start.pt == "br"){
			return(all.means)
			}		
		}	#end case for "tl" or "br"
	
	
	if(start.pt == "tl" || start.pt == "br"){

		diff.mat <- rotate.mat(mat1-mat2)
		# diff.mat[1:10, 1:10]
		
		min.idx <- min(diff.mat)
		max.idx <- max(diff.mat)
		all.idx <- min.idx:max.idx
		
		for(i in 1:length(all.idx)){
			report.progress(i, length(all.idx))
			idx.locale <- which(diff.mat == all.idx[i])
			diag.v <- mat[idx.locale]
			all.means[i] <- mean(diag.v, na.rm = TRUE)
			}

		
		if(start.pt == "tl"){
			return(all.means)
			}
		if(start.pt == "br"){
			return(rev(all.means))
			}		
		}
			
	} #end function