\name{as.matrix.Codelink}
\alias{as.matrix.Codelink}
\title{Return a matrix of intensity values}
\description{
	Takes a Codelink object and returns a matrix with
	the intensity values available.
}
\usage{
  \method{as.matrix}{Codelink}(x, ...)
}
\arguments{
  \item{x}{an object of class "Codelink".}
  \item{...}{additional arguments added to generic as.matrix since R-2.5.x}
}
\value{A matrix with the intensity values.}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{\code{\link{as.matrix}}}
\examples{
\dontrun{
	data(codelink.example)
	mat <- as.matrix(codelink.example)
	is(mat)
}
}
\keyword{documentation}
\keyword{utilities}
