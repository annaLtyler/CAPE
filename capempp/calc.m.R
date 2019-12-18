calc.m <- function(markers,beta.m,se,beta.cov) {
	# print(markers)		
	beta.main <- beta.m[,1:2]
	beta.inter <- beta.m[,3]
	n.rows <- dim(beta.main)[1] #No of ETs
	n.cols <- dim(beta.main)[2]

	se.main <- se[,1:2]
	se.inter <- se[,3]

	
	if(n.rows == n.cols){
		act_delta <- solve(beta.main)%*%beta.inter
		}else{
		act_delta <- try(pseudoinverse(beta.main)%*%beta.inter, silent = TRUE)
		if(class(act_delta) == "try-error"){
			act_delta <- c(NA, NA)
			}
		}
	
	m12 = act_delta[1]/(1+act_delta[2])
	m21 = act_delta[2]/(1+act_delta[1])
	

	results <- cbind(markers[1],markers[2],m12,m21)
	colnames(results) <- c("marker 1","marker 2","m12","m21")
	rownames(results)  <- NULL
	return(results)
	
	}