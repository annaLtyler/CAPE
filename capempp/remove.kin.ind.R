#this function deletes individuals from the kinship object
#so that it matches the data.obj

remove.kin.ind <- function(data.obj, kin.obj){
	
	if(class(kin.obj) == "matrix"){
		ind.locale <- match(rownames(data.obj$pheno), rownames(kin.obj))
		new.kin.obj <- kin.obj[ind.locale, ind.locale]
		}else{
		new.kin.obj <- lapply(kin.obj, function(x) 
			x[match(rownames(data.obj$pheno), rownames(x)), match(rownames(data.obj$pheno),colnames(x))])
		}
	
	return(new.kin.obj)	
	}