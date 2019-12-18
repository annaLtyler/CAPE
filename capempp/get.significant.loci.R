#This function takes in a cross object and returns a
#table of the significant loci for the 1D scan. The
#default behavior is to get the loci that clear the
#minimum threshold listed in the data objed, but this
#can be overridden by specifying a threshold.

get.significant.loci <- function(data.obj, threshold = NULL){
	
	cur.dir <- getwd()
		
	#This function translates the significant hits
	#into a more human readable format
	translate.sig.num <- function(sig.row){
		named.results <- NULL
		for(i in 1:length(sig.row)){
			named.results <- c(named.results, dimnames(stat.table)[[i]][sig.row[i]])
			}
		#add on the t statistic
		named.results <- c(named.results, round(stat.table[sig.row[1], sig.row[2], sig.row[3]], 3))
		#the chromosome number
		snp.locale <- which(data.obj$geno.names[[3]] == named.results[1])
		chrom <- data.obj$chromosome[snp.locale]
		named.results <- c(named.results, chrom)
		#and the chromosomal location
		loc <- data.obj$marker.location[snp.locale]
		named.results <- c(named.results, loc)
		return(named.results)
		}


	stat.table <- data.obj$oneDscan.t.stats	
	
	#find the minimum threshold value,
	#pull out the significant loci and
	#translate the table.
	if(is.null(threshold)){
		min.thresh <- get.threshold(data.obj, min)
		}else{
			min.thresh <- threshold #if a threshold is specified, use it instead
			}
	
	
	sig <- which(abs(stat.table) >= min.thresh, arr.ind = TRUE)
		
	
	if(length(sig) > 0){
		translated.table <- t(apply(sig, 1, translate.sig.num))
		#convert the X to a number, sort, and convert back
		x.locale <- which(translated.table[,5] == "X")
		translated.table[x.locale,5] <- 20
		sorted.table <- sort.by.then.by(translated.table, sort.cols = c(2,5), col.type = c("c", "n"), decreasing = FALSE)
		x.locale <- which(sorted.table[,5] == 20)
		sorted.table[x.locale,5] <- "X"
			
		}else{
			sorted.table <- matrix(rep("no.sig.loci", 6), nrow = 1)
			}
			
			
	colnames(sorted.table) <- c("locus", "phenotype", "allele", "t.stat", "chromosome", "position")
	setwd(cur.dir)
	
	return(sorted.table)

}
