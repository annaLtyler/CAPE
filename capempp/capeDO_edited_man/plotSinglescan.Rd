\name{plotSinglescan}
\alias{plotSinglescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Plot the results of singlescan
}
\description{
This function plots the results obtained from the single-marker regression performed by singlescan. The effects (\eqn{\beta}) of each regression on each phenotype or eigentrait are plotted as a vertical line. Chromosomes and traits for plotting can be specified.
}
\usage{
plotSinglescan(data.obj, chr = NULL, traits = NULL, 
standardized = TRUE, view = c("overview", "detailed"), 
allele.labels = NULL, include.covars = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{capeDO}}.
}
	\item{chr}{
		An optional vector indicating which chromosomes should be plotted. If NULL, the default, all chromosomes are plotted.
	}
	\item{traits}{
		An optional vector indicating which traits should be plotted. If NULL, the default, all traits are plotted.
	}
	\item{standardized}{
		A logical value. If TRUE, the absolute value of the regression t statistics (\eqn{\beta/\sigma}) are plotted. If FALSE, the raw regression coefficients (\eqn{\beta}) are plotted.
	}
	\item{view}{
		A character vector specifying whether the plot should be stretched to show a detailed view of the allele effects on each chromosome.
	}
	\item{allele.labels}{
		An optional vector specifying the name of each allele for the legend.	
		}
	\item{include.covars}{
		A logical value specifying whether the covariates should be included in the plot. In some cases, covariates can have a much stronger effect than any of the genetic markers and make it difficult to see the relative effects of the markers.
	}

}
\value{
Nothing is returned from this function. It produces a plot of the allele effects. If the effects are not standardized, the plot is split into upper and lower panels. The upper panel shows the LOD score for each locus overall. The lower panel shows the individual effect size of each allele. If the effects are standardized by dividing by the standard error, there is only one panel in the plot and the standardized effect of each allele is shown. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{singlescan}}, \code{\link{select.markers.for.pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}

