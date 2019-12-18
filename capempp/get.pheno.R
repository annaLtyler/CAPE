get.pheno <- function(data.obj, scan.what = c("eigentraits", "normalized.traits", "raw.traits"), covar = NULL){
	
	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentraits, otherwise, use raw phenotypes
	scan.what <- scan.what[1]
	is.ET <- c(grep("eig", scan.what, ignore.case = TRUE), grep("ET", scan.what, ignore.case = TRUE))
	is.raw <- grep("w", scan.what, ignore.case = TRUE)	 
	is.norm <- grep("o", scan.what, ignore.case = TRUE)	 
	

	if(length(is.ET) > 0){ 
		el.idx <- which(names(data.obj) == "ET")
		if(length(el.idx) == 0){
			stop("There are no eigentraits. Run get.eigentraits() to generate eigentraits.")
			}
		}
		
	if(length(is.raw) > 0){
		el.idx <- which(names(data.obj) == "pheno")
		}

	if(length(is.norm) > 0){
		el.idx <- which(names(data.obj) == "raw.pheno")
		if(length(el.idx) == 0){
			el.idx <- which(names(data.obj) == "pheno")
			}
		}
	
	pheno <- data.obj[[el.idx]]
	
	if(!is.null(covar)){
		covar.info <- get.covar(data.obj)
		covar.locale <- which(covar.info$covar.names %in% covar)
		models <- apply(pheno, 2, function(x) lm(x~covar.info$covar.table[,covar.locale,drop=FALSE]))
		resids <- lapply(models, residuals)
		resid.table <- t(list2Matrix(resids))
		return(resid.table)
		}
	
	return(pheno)
	
	
}