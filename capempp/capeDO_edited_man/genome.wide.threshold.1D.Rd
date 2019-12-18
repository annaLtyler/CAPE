\name{genome.wide.threshold.1D}
\alias{genome.wide.threshold.1D}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Calculate a genome-wide significance threshold for the single-variant scan
}
\description{
This function uses permtuation testing to calculate a genome-wide significance threshold for the single-variant scan. Two user-defined thresholds are calculated: the lower threshold (higher alpha) can be used to determine which variants will be used in the pairwise scan. The higher threshold (lower alpha) is used to determine which variants are used as covariates in the pairwise scan. In each permutation, the phenotype or eigentrait is shuffled and all markers are retested with the permuted phenotype. The regression coefficients are collected from each permutation and the extreme value distribution is used to determine thresholds for the user-defined alpha values.
}
\usage{
genome.wide.threshold.1D(data.obj, n.perm = 100, 
scan.what = c("eigentraits", "raw.traits"), ref.allele = NULL, 
alpha.for.pairs = 0.05, alpha.for.covar = 0.01, verbose = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{capeDO}}.
	}
  \item{n.perm}{
	The number of permutations to perform. The default is 100.
	}
	\item{scan.what}{
	A character value that uniquely specifies whether the eigentraits or phenotypes should be scanned.
	}
	\item{ref.allele}{
		\code{\link{capeDO}} requires that one of the alleles in the population be selected as a reference allele. The effects of all alleles are then reported as the effects relative to the reference allele. In the DO population B represents the Black6 mouse. Because this strain is considered a very standard strain, the default reference allele is B. If another allele has a positive effect, it means that it increases the given phenotype relative to that in the Black6 mouse. Any other allele can be selected as the reference allele simply by specifying it with this argument. 
	}
  \item{alpha.for.pairs}{
	The alpha value used to calculate the threshold for which variants will be used in the pair scan. In \code{pairscan} the use of this threshold is optional.
}
  \item{alpha.for.covar}{
	The alpha value used to calculate the threshold for which variants will be used as covariates in the pairwise scan.
}
  \item{verbose}{
	A logical value specifying whether the progress of the scan should be printed to the screen. A value of TRUE causes the progress to be printed to the screen.
}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{singlescan}}
\code{\link{pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

