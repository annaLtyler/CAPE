\name{bin.vector}
\alias{bin.vector}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Bins a vector into values near the specified values.
}
\description{
This function converts a vector of continuous values into a vector of specified values by moving each continuous vector to the nearest value in the vector specified by the user. For example, a vector of values ranging from 0 to 1 can be binned into a vector of only 0, 0.5 and 1. Each value in the continuous vector is rounded to the nearest of the desired values. This is useful in converting continuous phenotypes into the traditional 0, 0.5 and 1 values.
}
\usage{
bin.vector(vectorX, bins = seq(0,1,0.5))
}

\arguments{
  \item{vectorX}{
	The vector to be binned.
	}
	\item{bins}{
		A vector of the only values desired in the final vector.
	}
}
\value{
A vector of the same length as the original vector is returned. Each value has been rounded to the nearest value specified by the argument bins.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{arith}% __ONLY ONE__ keyword per line
\keyword{internal}