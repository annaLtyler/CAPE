#This function plots the 36-state boxplot of phenotypes
#given a marker name.


pheno.boxes <- function(data.obj, geno.obj = NULL, marker, pheno = NULL, sort.by = c("genotype", "median", "mean")){
	
	sort.check <- grep("g", sort.by)
	if(length(sort.check) > 0){
		sort.by = "genotype"
		}
	
	
	geno <- get.geno(data.obj, geno.obj)
	marker.locale <- which(data.obj$geno.names[[3]] == marker)
	if(length(marker.locale) == 0){
		stop(paste("I couldn't find marker:", marker))
		}
	
	marker.geno <- geno[,,marker.locale]
	binned.geno <- apply(marker.geno, 2, bin.vector)
	
	get.state <- function(one.geno){
		alleles.which <- names(which(one.geno != 0))
		if(length(alleles.which) == 1){
			return(paste(alleles.which, alleles.which, sep = "/", collapse = ""))
			}
		if(length(alleles.which) == 2){
			sorted.alleles <- sort(alleles.which)
			return(paste(sorted.alleles[1], sorted.alleles[2], sep = "/"))
			}
		}
	
	all.alleles <- dimnames(geno)[[2]]
	all.states <- apply(pair.matrix(all.alleles, self.pairs = TRUE), 1, function(x) paste(x[1], x[2], sep = "/"))
	sorted.states <- sort(all.states)
	
	#get the genotype state of each animal
	ind.states <- apply(binned.geno, 1, get.state)
	
	if(is.null(pheno)){
		pheno.locale <- 1:dim(data.obj$pheno)[2]
		pheno.names <- colnames(data.obj$pheno)
		}else{
		pheno.locale <- get.col.num(data.obj$pheno, pheno)
		pheno.names <- pheno
		}

	layout.mat <- matrix(1:length(pheno.locale), ncol = 1)
	layout(layout.mat)
	for(ph in 1:length(pheno.locale)){
		pheno.vals <- data.obj$pheno[,pheno.locale[ph]]
		
		pheno.list <- vector(mode = "list", length = length(sorted.states))
		names(pheno.list) <- sorted.states
		for(i in 1:length(pheno.list)){
			val.locale <- which(ind.states == sorted.states[i])
			pheno.list[[i]] <- pheno.vals[val.locale]
			}
		
		if(sort.by == "mean"){
			list.order <- order(unlist(lapply(pheno.list, function(x) mean(x, na.rm = TRUE))))
			sorted.list <- pheno.list[list.order]
			}
		if(sort.by == "median"){
			list.order <- order(unlist(lapply(pheno.list, function(x) median(x, na.rm = TRUE))))
			sorted.list <- pheno.list[list.order]			
			}
		if(sort.by == "genotype"){
			sorted.list <- pheno.list
			}	

		boxplot(sorted.list, main = paste(pheno[ph], "\n", marker), names = NA)
		plot.dim <- par("usr")
		yrange = plot.dim[4]-plot.dim[3]
		par(xpd = TRUE)
		text(x = 1:length(sorted.list), y = plot.dim[3]-(yrange*0.06), names(sorted.list), srt = 90)
		par(xpd = FALSE)
		}
	
	
}

