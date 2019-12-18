make.fake.kin <- function(data.obj, geno.obj, type = c("overall", "ltco")){
	
	geno <- get.geno(data.obj, geno.obj)
	kin <- cov(t(geno[,2,]))
	kin.list <- kin + abs(min(kin))
	
	if(type == "ltco"){
		chr <- unique(data.obj$chromosome)
		chr.pairs <- pair.matrix(chr, self.pairs = TRUE)
		chr.names <- apply(chr.pairs, 1, function(x) paste(x, collapse = ","))
		kin.list <- lapply(1:nrow(chr.pairs), function(x) kin)
		names(kin.list) <- chr.names	
		}
	
	return(kin)
	
	}