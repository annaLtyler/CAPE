#This function plots phenotypes in correlation with a supplied vector
#for example, percentage strain background or first PC of genotype matrix


plot.pheno.V <- function(data.obj, v, v.label, covar, pdf.label = "Correlation.of.Phenotypes.with.Vector.pdf"){
	
	pheno <- get.pheno(cross, "raw.traits", covar = covar)
	layout.mat <- get.layout.mat(ncol(pheno))
	pdf(pdf.label, width = ncol(layout.mat)*3, height = nrow(layout.mat)*3)
layout(layout.mat)
for(i in 1:ncol(pheno)){
	model <- lm(pheno[,i]~v)
	pval <- signif(coefficients(summary(model))[2,4], 2)
	plot(v, pheno[,i], main = paste(colnames(pheno)[i], "\np =", pval), xlab = "Percent DBA", ylab = colnames(pheno)[i])	
	abline(model)
	}
	dev.off()

	
}