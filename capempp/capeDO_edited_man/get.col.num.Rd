\name{get.col.num}
\alias{get.col.num}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Use column titles to retrieve column numbers
}
\description{
This function is used internally to tranlate column names to column numbers.
}
\usage{
get.col.num(data.mat, col.which = NULL, dim.which = 2)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.mat}{
	A matrix with column names
}
  \item{col.which}{
	The names or numbers of the columns for which numbers are wanted.
}
  \item{dim}{
	A numeric value indicating the dimension of the array in which to look for labels.
}
}
\examples{

mat <- matrix(rnorm(50), ncol = 5)
colnames(mat) <- letters[1:5]
print(mat)
col.num <- get.col.num(mat, c("c", "e"))
print(col.num)

col.num <- get.col.num(mat, c("b", "g"))
print(col.num)

}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

