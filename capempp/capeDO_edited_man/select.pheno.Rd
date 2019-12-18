\name{select.pheno}
\alias{select.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Select phenotypes for analysis
}
\description{
This function is used to pare down a matrix of phenotypes to only those that will be included in the analysis.
}
\usage{
select.pheno(data.obj, pheno.which, min.entries = 5, 
scale.pheno = FALSE, rank.norm.pheno = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{capeDO}}.
	}
  \item{pheno.which}{
	A vector of either column numbers of the desired phenotypes, or character strings that uniquely identify the columns containing the desired phenotypes.
	}
	\item{min.entries}{
		The miniminum number of values allowed for a phenotype. Phenotypes with fewer than this number of values will not be selected for analysis.
	}
	\item{scale.pheno}{
		A logical value indicating whether the phenotypes should be mean centered and standardized.
	}
	\item{rank.norm.pheno}{
		A logical value indicating whether the phenotypes should be normalized using rank Z normalization.
	}
}
\value{
This function returns the data object with a phenotype matrix that contains only the specified phenotypes.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}

