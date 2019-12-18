#This function bins the entire genome into linkage 
#blocks based on pairwise genotype correlation
#first it calculated pairwise correlations between
#markers on a single chromosome by sliding a window
#of a given depth across the chromosome. The density
#of correlations calculated can be adjusted using min.cor
#and overlap. min.cor is for determining the window width
#it is based on how far away markers get from a reference
#before their correlation is the minimum correlation.
#overlap gives the percentage of overlap between windows.
#higher overlaps give smoother sampling.

# setwd("~/Documents/Data/DO/Results/all_pheno")
# source('~/Documents/git_repositories/useful_r_code/pair.matrix.R')
# source('~/Documents/git_repositories/capeDO/consec.pairs.R')
# data.obj <- readRDS("crossDO.RData")
# geno.obj <- readRDS("~/Documents/Data/DO/data_for_cape/capeDOGeno2016-11-22.RData")

# setwd("~/Documents/Data/DO/Results/Genome_Binning")

full.geno.linkage.blocks <- function(data.obj, geno.obj, pair.depth = 400, binning.window.size = 20, plot.blocks = TRUE, save.calculations = TRUE){
	

require(igraph)
require(pbapply)


#================================================================
#internal functions
#================================================================
get.cor <- function(pair.num, chr.geno){
	marker1.geno <- as.vector(chr.geno[,,pairs.to.test[pair.num,1]])
	marker2.geno <- as.vector(chr.geno[,,pairs.to.test[pair.num,2]])
	return(cor(marker1.geno, marker2.geno, use = "complete"))
	}
	
#================================================================

	
	all.ch <- data.obj$chromosome
	u_chr <- unique(all.ch)
	
	for(ch in 1:length(u_chr)){
	
		
		cat("Chromosome", u_chr[ch], "\n")		

		cat("\tFinding marker pairs for correlation calculations...\n")
		all.idx <- which(data.obj$chromosome == u_chr[ch])
		chr.geno <- geno.obj[,,all.idx]
		
		pairs.to.test <- pair.idx(num = length(all.idx), depth = pair.depth)

		num.pairs <- nrow(pairs.to.test)
		cat("\tCalculating correlations for", num.pairs, "pairs\n")
				
 		if(!file.exists(paste0("adj.mat.Chr", u_chr[ch], ".RData"))){
			cat("\tCalculating pairwise correlations...\n")
			all.block.cor <- pbsapply(1:num.pairs, function(x) get.cor(x, chr.geno))
 	
 			net <- graph.edgelist(pairs.to.test, directed = FALSE)
			E(net)$weight = all.block.cor
				
			cat("\tLinearizing pairwise correlations...\n")
			adj.mat <- as.matrix(as_adjacency_matrix(net, type = "both", attr = "weight"))
			adj.mat[which(adj.mat == 0)] <- NA
		
			if(save.calculations){
				saveRDS(adj.mat, paste0("adj.mat.Chr", u_chr[ch], ".RData"))
				}
			}else{
			adj.mat <- readRDS(paste0("adj.mat.Chr", u_chr[ch], ".RData"))
			}
		
		if(!file.exists(paste0("Diag.Means.Chr", u_chr[ch], ".RData"))){
			diag.means <- diagonal.means(adj.mat, "tl")
			
			if(save.calculations){
				saveRDS(diag.means , paste0("Diag.Means.Chr", u_chr[ch], ".RData"))
				}
			}else{
				diag.means <- readRDS(paste0("Diag.Means.Chr", u_chr[ch], ".RData"))
				}
		
		if(plot.blocks){
			pdf(paste0("Peaks.Binned.Chr",u_chr[ch], ".pdf"))
			smoothed.peaks <- bin.curve(diag.means, plot.peaks = TRUE, window.size = binning.window.size, amp.min = 0.01)
			peaks <- smoothed.peaks$bins
			dev.off()
			}
		
		
		# start.time <- Sys.time()		
		if(plot.blocks){
			color.breaks <- c(0.25, 0.5, 0.75)
			cat("\tPlotting Linkage blocks...\n")
			jpeg(paste0("Linkage.Blocks.Chr.", u_chr[ch], ".jpg"), 1000, 1000)
			plot.adj.mat.with.blocks(adj.mat, l.blocks = peaks, plot.label = paste("Linkage Blocks Chromosome", u_chr[ch]), color.breaks = color.breaks)
			dev.off()
			}



		cat("\tAssigning markers to linkage blocks...\n")
		
		#assign the markers to linkage blocks based on the peaks
		chr.blocks <- assign.linkage.blocks(data.obj, u_chr[ch], peaks)		
		
		# block.sizes <- get.block.sizes(chr.blocks)
		if(ch == 1){
			all.linkage.blocks <- chr.blocks
			names(all.linkage.blocks) <- paste0("Chr1_", 1:length(all.linkage.blocks))
			}else{
			start.pos <- length(all.linkage.blocks)+1
			end.pos <- start.pos + length(chr.blocks) - 1
			names(chr.blocks) <- paste0("Chr", u_chr[ch], "_", 1:length(chr.blocks))
			all.linkage.blocks[start.pos:end.pos] <- chr.blocks
			}
		

	} #end looping through chromosomes
	
	return(all.linkage.blocks)
	
}
	
