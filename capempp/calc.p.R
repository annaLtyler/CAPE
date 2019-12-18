calc.p <- function(data.obj, pval.correction = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")) {
	
	
	pval.correction = pval.correction[1]
	data.obj$pval.correction <- pval.correction	
	
	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm
		
	if(is.null(influences.org)){
		stop("error.prop() with perm = FALSE must be run before running calc.p()")
		}
		
	if(is.null(influences.perm)){
		warning("error.prop() with perm = TRUE must be run before running calc.p()")
		}

	n.gene <- dim(data.obj$geno.for.pairscan)[2] #get the number of genes used in the pair scan
	n.pairs <- dim(influences.org)[1] #the number of pairs scanned in the pairscan
        	
    marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
    colnames(marker.mat) <- c("marker1", "marker2")


	#### Combinine across permutations#####
	#get the t statistics for all permutations
	m12.null <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	m21.null <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	m.null <- c(m12.null, m21.null)

	m12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
	m21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])
	m <- c(m12, m21)

	#changed calculation of p value to account for the asymmetric m12/m21 distribution
	#I now calculate the p value based on above and below 0 separately
	low.null <- m.null[which(m.null < 0)]; low.fun <- ecdf(abs(low.null))
	high.null <- m.null[which(m.null > 0)]; high.fun <- ecdf(high.null)

	all.emp.p <- matrix(NA, ncol = 1, nrow = length(m))
	high.m.locale <- which(m > 0)
	low.m.locale <- which(m < 0)

	low.p <- 1-low.fun(abs(m[low.m.locale]))
	high.p <- 1-high.fun(abs(m[high.m.locale]))

	all.emp.p[low.m.locale,1] <- low.p
	all.emp.p[high.m.locale,1] <- high.p
	# hist(all.emp.p)

	m12 <- matrix(c(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4]))), ncol = 5)	
	colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE")	

	m21 <- matrix(c(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6]))), ncol = 5)
	colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE")

	p.adjusted <- p.adjust(all.emp.p, method = pval.correction)
	final.table <- rbind(m12, m21)
	final.table <- cbind(final.table, all.emp.p, p.adjusted)
	colnames(final.table)[c(6,7)] <- c("P_empirical", "p.adjusted")
	final.table <- final.table[order(as.numeric(final.table[,"|Effect|/SE"]), decreasing = TRUE),]

	# head(final.table)

	data.obj$var.to.var.p.val <- final.table
	
	return(data.obj)
}
