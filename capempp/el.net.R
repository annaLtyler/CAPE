el.net <- function(x, y, num.alpha = 10, alpha = NULL, stand = TRUE, plot.diagnostics = FALSE, pdf.name = NULL){
	
	if(!is.null(alpha) && num.alpha > 1){
		stop("If alpha is specified, num.alpha must be 1.")
		}
	
	if(is.null(pdf.name)){pdf.name = "el.net.models.pdf"}
	
	# Find optimal alpha
	# a = seq(1/num.alpha, 1, 1/num.alpha)
	if(num.alpha == 1){
		a <- alpha
		}else{
		a = 2^seq(from = -(num.alpha-1), to = 0, 1)
		}

	cvm.mat <- matrix(NA, nrow = num.alpha, ncol = 100)
	lambda.mat <- matrix(NA, nrow = num.alpha, ncol = 100)

	opt.cvm = Inf
	opt.idx = NaN
	if(plot.diagnostics){pdf(pdf.name)}
	for(i in 1:length(a)){
		curr.mod = cv.glmnet(x, y, intercept = FALSE, alpha = a[i], standardize = stand)
		cvm.mat[i,1:length(curr.mod$cvm)] <- curr.mod$cvm
		lambda.mat[i,1:length(curr.mod$lambda)] <- curr.mod$lambda
		
		if(plot.diagnostics){plot(curr.mod); mtext(text = paste("alpha =", a[i], ":: min cross-validation error =", signif(min(curr.mod$cvm), 2)), outer = TRUE, line = -1.2)}
		# print(min(curr.mod$cvm))
		# print(min(curr.mod$cvm) < opt.cvm)
		
		if(min(curr.mod$cvm) < opt.cvm){
			opt.cvm = min(curr.mod$cvm)
			opt.idx = i
		}
	}
	if(plot.diagnostics){dev.off()}
	# Find optimal lambda

	if(plot.diagnostics){
		#first we need to interpolate the cvm.mat at fixed values of lambda
		min.lambda <- min(log10(lambda.mat), na.rm = TRUE)
		max.lambda <- max(log10(lambda.mat), na.rm = TRUE)
		lambda.seq <- seq(min.lambda, max.lambda, length.out = 1000)
		interp.cvm <- matrix(NA, nrow = num.alpha, ncol = 1000)
		for(i in 1:num.alpha){
			interp.fun <- approxfun(log10(lambda.mat[i,]), cvm.mat[i,], method = "constant", rule = 1)
			interp.cvm[i,] <- interp.fun(lambda.seq)
			}
		pdf(paste("cvm.mat.", pdf.name, sep = ""))
		image(rotate.mat(interp.cvm))
		dev.off()
		}

	curr.mod = cv.glmnet(x, y, intercept = FALSE, alpha = a[opt.idx], standardize = stand)
	
	opt.mod = glmnet(x, y, intercept = FALSE, alpha = a[opt.idx], standardize = stand, lambda = curr.mod$lambda.min)
	
	return(list(a[opt.idx], opt.mod, curr.mod$cvsd))
	
}