\name{capeDO-package}
\alias{capeDO-package}
\alias{capeDO}
\docType{package}
\title{ 
	Combinatorial Analysis of Epistasis and Pleiotropy for the Diversity Outbred Mice
}
\description{
This package infers predictive networks between genetic variants and between genetic variants and phenotypes. It is specifically designed for use with the multi-allele mouse population known as the Diversity Outbred (DO) mice. It uses complementary information of pleiotropic gene variants across different phenotypes to resolve models of epistatic interactions between alleles. To do this, cape reparameterizes main effect and interaction coefficients from a pairwise variant regressions into directed influence parameters. These parameters describe how alleles each other, in terms of suppression and enhancement, as well as how gene variants influence phenotypes. In this framework, alleles originating from individual founders are selected for pairwise regressions based on their main effects on individual phenotypes. These alleles are then treated as independent genetic markers in the pairwise scans. The result is different founder alleles from the same locus may interact with each other in the form of intra-locus epistasis, or dominance. All of the final interactions are reported in terms of founder alleles at a given genomic location.
}
\details{
\tabular{ll}{
Package: \tab cape\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-12-30\cr
License: \tab GPL-3\cr
}
The cape analysis begins by reading in a genetic data set with three functions: \code{\link{read.pheno}}, \code{\link{read.geno}}, and \code{\link{combine.geno.pheno}}. These functions build a data object including a three-dimensional genotype array, phenotype data, and information about the genetic markers. Throughout the analysis, results are added to the data object. This data object is an argument in most functions, and is referred to as \code{data.obj}. Because this package uses pleiotropy to resolve models of epistasis, the phenotypes used should have common underlying molecular players, but not be perfectly correlated. In general phenotypes should be correlated with a Pearson r between 0.4 and 0.8. The phenotypes of interest are then decomposed into eigentraits using \code{\link{get.eigentraits}}. Any number of phenotypes can be decomposed, but \code{cape} requires between two and 12 eigentraits for the analysis. The phenotype decomposition into eigentraits maximizes the complementary information between the phenotypes. Before investigating epistatic interactions through a pair-wise scan of the genetic variants, a single-variant scan (\code{\link{singlescan}}) is run. This allows for thresholding of markers for the pair scan based on main effect size. The single-variant scan also allows selection of variants with very large main effects to be used as covariates in the pair scan. The pair scan (\code{\link{pairscan}}) performs a regression on each variant pair. Finally, the coefficients from the pairwise scan are reparametrized to yield directional influences between variants and from variants to phenotypes. In this DO version of \code{cape}, the variants contain founder allele information in addition to genomic location information. The final result is an asymmetric adjacency matrix describing how variants from the different founders influence each other to affect phenotype. The p values of these influences are corrected for multiple tests.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
Maintainer: Anna L. Tyler <anna.tyler@jax.org>
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010.

Churchill, G. A., Gatti, D. M., Munger, S. C., & Svenson, K. L. (2012). The diversity outbred mouse population. Mammalian Genome, 23(9-10), 713-718.
}

\keyword{package}

