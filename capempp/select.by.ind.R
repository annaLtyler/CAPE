#This script subsets a cross by a given 
#column. It takes in a cross object and
#an expression indicating the subset of
#interest. This is in the following form:
#column.name == value, or column.name != value
#or calumn.name < value, etc.
#where colum.name is an unambiguous string
#identifying a particular column
#The script returns the subsetted cross.

select.by.ind <- function(data.obj, geno.obj = NULL, geno.or.pheno = "pheno", expr){
	
	geno <- get.geno(data.obj, geno.obj)
	
	mat.to.sub.test <- grep("p", geno.or.pheno)
	if(length(mat.to.sub.test) > 0){
		sub.mat <- data.obj$pheno
		}else{
			sub.mat <- geno
			}
	
	expr.pieces <- strsplit(expr, "\\ ")
	if(length(expr.pieces[[1]]) != 3){
		stop("Expression must be in the format 'colname comparison value'")
		}
	
	if(geno.or.pheno == "geno"){
		col.locale <- which(data.obj$geno.names[[3]] == as.character(expr.pieces[[1]][1]))
		}else{
		col.locale <- which(colnames(sub.mat) == as.character(expr.pieces[[1]][1]))
		}
	if(length(col.locale) == 0){
		stop("I can't find the column name: ", expr.pieces[[1]][1])
		}
		
	if(length(col.locale) > 1){
		stop("There is more than one column that matches the string: ", expr.pieces[[1]][1])
		}

	#if the final element in the expression is a number
	if(length(mat.to.sub.test) > 0){
		if(!is.na(as.numeric(expr.pieces[[1]][3]))){
			cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), as.numeric(expr.pieces[[1]][3]))
			}else{ #otherwise
				cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,col.locale]), expr.pieces[[1]][3])
				}
			}else{
			if(!is.na(as.numeric(expr.pieces[[1]][3]))){
			cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,1,col.locale]), as.numeric(expr.pieces[[1]][3]))
			}else{ #otherwise
				cl <- call(expr.pieces[[1]][2], as.numeric(sub.mat[,1,col.locale]), expr.pieces[[1]][3])
				}			
			}
			
	vals.locale <- which(eval(cl))	
	
	if(length(vals.locale) == 0){
		stop("There are no individuals that match this expression.")
		}
	
	if(length(vals.locale) < dim(data.obj$pheno)[1]){
		cat(dim(data.obj$pheno)[1] - length(vals.locale), "individuals were removed.\n")
		}
	
	data.obj$pheno <- data.obj$pheno[vals.locale,]
	data.obj$geno <- geno[vals.locale,,]
	if(!is.null(data.obj$raw.pheno)){
		data.obj$raw.pheno <- data.obj$raw.pheno[vals.locale,]
		}
	if(!is.null(data.obj$ET)){
		data.obj$ET <- data.obj$ET[vals.locale,]
		}	
	
	return(data.obj)
	
}