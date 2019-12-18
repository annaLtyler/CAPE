\name{read.geno}
\alias{read.geno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Read in and format genotype data for analysis by \code{\link{capeDO}}
}
\description{
This function reads in genotype data for \code{capeDO} analysis and formats it into a three-dimensional array. The array contains individuals in the first dimension, alleles in the second dimension and loci in the third dimension. This genotoype array is combined with the phenotype matrix created by \code{\link{read.pheno}}. The genotype and phenotype objects are combined using \code{\link{combine.geno.pheno}}.
}
\usage{
read.geno(path = ".", marker.info)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{path}{
	The path of the directory containing the genotype data. The default path is the current working directory.
	}
	\item{marker.info}{
	A matrix containing information about the genetic markers. The matrix must contain 3 columns. Column 1 contains the name of each marker, column 2 contains the chromosome number of each marker, and column 3 contains the chromosomal position of each marker.
	}
}
\details{
The genotype data for a multi-allele population are stored in a set of files rather than a single file that can contain the genotype data of a bi-allelic population. Only genotype files of the correct format should be contained in the genotype data directory. Each file should be named with the individual subject name. This function has been developed to work with the output of DOqtl, which estimates diplotype genotype probabilities for each founder allele in the diversity outbred population. Each file should contain a header with the diplotype genotype labeling each column. For example, "AA", "AB", "AC", etc.
Each row contains XXX
}
\value{
The genotype files are converted to a three-dimensional array that contains founder allele probabilities for each individual in the population. This array is combined with locus information from the marker.info matrix to create a genotype object containing three elements:
\item{geno}{A three-dimensional array containing founder allele probabilities for each individual at each locus.}
\item{chromosome}{A vector containing the chromosome on which each marker is found.}
\item{marker.location}{A vector containing the chromosomal position of each marker.}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{IO}

