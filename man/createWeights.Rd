\name{createWeights}
\alias{createWeights}
\alias{createWeights,CodelinkSet-method}
\alias{createWeights,Codelink-method}
\title{Create weight for Codelink chips}
\description{
  Create a weight matrix based on probe type and flags.
}
\usage{
	createWeights(object, type.weights = NULL, flag.weights = NULL)
}
\arguments{
	\item{object}{an object of class "Codelink" or "CodelinkSet".}
	\item{type.weights}{named character vector of type weights used.}
	\item{flag.weights}{named character vector of flag weights used.}
}
\author{Diego Diez}
\examples{
    data(codset)
    w <- createWeights(codset)
}
\keyword{documentation}
\keyword{utilities}
