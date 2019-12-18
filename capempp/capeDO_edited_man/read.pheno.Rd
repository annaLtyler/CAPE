\name{read.pheno}
\alias{read.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Read in and format the phenotype data for the analysis by \code{\link{capeDO}}
}
\description{
This function reads in phenotype data for \code{capeDO} analysis and formats it into an object used by \code{\link{combine.geno.pheno}}. A single comma-separated file containing only phenotype data is required.
Phenotype data must be numeric and contained in columns, with each individual occupying a row. A header 
indicating phenotype names is required.
}
\usage{
read.pheno(path = ".", filename = NULL, 
pheno.col = NULL, delim = ",", na.strings = "-")
}
%- maybe also usage for other objects documented here.
\arguments{
	\item{path}{
		The path of the directory containing the phenotype data. The default path is the current working directory.
		}
  \item{filename}{
	An optional character string with path name specifying the file to be read in. Omission of this argument will prompt a dialog box for selecting a file.
	}
  \item{pheno.col}{
	An optional numeric vector specifying which columns the phenotypes of interest are in. If omitted, all phenotypes are read in.
}
	\item{delim}{
		A character string indicating the delimeter in the data file. The default indicates a comma-separated file (",").
	}
	\item{na.strings}{
	The symbol used to denote missing data in the file. 	
	}
}
\details{
All phenotype data must be contained in a single file. The phenotypes should be listed in columns. Each row of the file corresponds to one individual. The file must contain a header indicating the names of the phenotypes.
}
\value{
The file is converted to a matrix that can be combined with a genotype object using \code{\link{combine.geno.pheno}}. The genotype object is read in with the function \code{\link{read.geno}}.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{IO}

