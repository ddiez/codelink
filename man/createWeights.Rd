\name{createWeights}
\alias{createWeights}
\title{Create weight for Codelink chips}
\description{
  Create a weight matrix to be used in limma.
}
\usage{
	createWeights(object, type=NULL)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{type}{weight assigned to each Probe\_type.}
}
\author{Diego Diez}
\examples{
\dontrun{
    data(codelink.example)
    w <- createWeights(codelink.example, type = list(FIDUCUAL = 0.1))
}
}
\keyword{documentation}
\keyword{utilities}
