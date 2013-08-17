\name{createWeights}
\alias{createWeights}
\title{Create weight for Codelink chips}
\description{
  Create a weight matrix based on probe type and flags.
}
\usage{
	createWeights(object, type.weights = NULL, flag.weights = NULL)
}
\arguments{
	\item{object}{an object of class "Codelink" or "CodelinkSet".}
	\item{type.weights}{type weights used.}
	\item{flag.weights}{flag weights used.}
}
\author{Diego Diez}
\examples{
\dontrun{
    data(codelink.example)
    w <- createWeights(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
