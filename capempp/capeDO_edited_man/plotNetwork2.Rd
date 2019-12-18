\name{plotNetwork2}
\alias{plotNetwork2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the final epistatic network in a traditional network view.
}
\description{
This function plots the final results in a layout different to both \code{\link{plotVariantInfluences}} and \code{\link{plotNetwork}}. In this view, the network is plotted with a traditional network layout. The positional information in \code{\link{plotNetwork}} is lost, but in this view it is easier to see the structure of the overall network in terms of hub and peripheral nodes. In this view, each node is plotted as a pie-chart, and the main effects of the node are indicated as positive (green), negative (red), or not-significant (gray). Significant interactions are shown as green (positive) or red (negative) arrows between nodes. The function \code{\link{get.network}} must be run before plotting the network.
}
\usage{
plotNetwork2(data.obj, p.or.q = 0.05,  collapsed.net = FALSE, 
r2.thresh = 0.5, layout.matrix = NULL, node.radius = 1, 
label.nodes = TRUE, label.offset = 0, label.cex = 1, 
legend.radius = 1, legend.cex = 1, arrow.offset = 1, 
edge.lwd = 2)
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
\code{\link{get.network}}, \code{\link{plotVariantInfluences}}, \code{\link{plotNetwork}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

