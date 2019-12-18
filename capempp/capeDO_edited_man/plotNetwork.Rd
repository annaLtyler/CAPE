\name{plotNetwork}
\alias{plotNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the final epistatic network
}
\description{
This function plots the final results with a using a different layout than \code{\link{plotVariantInfluences}}. Instead of the adjacency matrix, this function plots interactions between markers with arrows indicating the direction of influence. The function \code{\link{get.network}} must be run before plotting the network.
}
\usage{
plotNetwork(data.obj, trait = NULL, chr = NULL, 
show.effect.size = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{capeDO}}.
	}
  \item{trait}{
	A character vector indicating which traits should be plotted. If NULL, all traits are plotted.
}
  \item{chr}{
	A character vector indicating which chromosomes should be plotted. If NULL, all chromosomes are plotted.
	}
\item{show.effect.size}{
	A logical value indicating whether the height of chromosome main effect bars should vary with the effect size.
	}
}
\seealso{
\code{\link{get.network}}, \code{\link{plotVariantInfluences}}, \code{\link{plotNetwork2}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

