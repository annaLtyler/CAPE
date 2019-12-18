# pop <- readRDS("~/Documents/Grants/R21/Scleroderma/Results/capeRel_testing/pop.RData")
# geno <- readRDS("~/Documents/Grants/R21/Scleroderma/Results/capeRel_testing/popGeno.RData")
# lig4.snps <- get.snps.in.gene("LIG4", "human", 5e5, 5e5)
# snp.locale <- which(geno$marker.names %in% lig4.snps[,1])



get.snps.near.snp <- function(snp.name, organism = c("mouse", "human"), upstream.buffer = 0, downstream.buffer = 0){
	require(biomaRt)
		if(organism == "mouse"){		
		# lib <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
			lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2017.archive.ensembl.org")
			snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="mmusculus_snp", host = "may2017.archive.ensembl.org")	
		}else{
		# lib <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")	
		lib <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2017.archive.ensembl.org")
		snp.db = useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = "may2017.archive.ensembl.org")
		}

	snp.coord <- getBM(c("chr_name", "chrom_start", "chrom_end"), "snp_filter", values = snp.name, mart = snp.db)
	nearby.snps <- getBM("refsnp_id", "chromosomal_region", paste(snp.coord[1,1], as.numeric(snp.coord[1,2] - upstream.buffer), as.numeric(snp.coord[1,3]+downstream.buffer), sep = ":"), snp.db)
	return(nearby.snps)
	
	
}