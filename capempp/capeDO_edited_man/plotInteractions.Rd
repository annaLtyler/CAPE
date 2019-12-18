\name{plotInteractions}
\alias{plotInteractions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate interaction plots for the significant interactions.
}
\description{
This function plots interaction plots for all significant interactions and all phenotypes and eigentraits.
}
\usage{
plotInteractions(data.obj, p.or.q = 0.05,
pdf.label = "Interaction.Plots.pdf", raw.pheno = FALSE,
error.bars = TRUE, covar = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{capeDO}}.
}
	\item{p.or.q}{
		A threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant. Only marker pairs with p or q values below this threshold will be plotted.
	}
	\item{pdf.label}{
		A character string indicating the name of the file to which the interactions are plotted.
	}
	\item{raw.pheno}{
		A logical value indicating whether the raw phenotypes or normalized phenotypes should be used in the plots.
	}
	\item{error.bars}{
		A logical value indicating whether error bars (here standard error) should be plotted on the interaction plots.
	}
	\item{covar}{
		An optional vector indicating which which variables should be used as covariates to correct the resulting plots.
		}

}
\value{
No values are returned. The function generates a interaction plots that are saved to a pdf.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

