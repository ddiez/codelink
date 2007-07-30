\name{dim.Codelink}
\alias{dim.Codelink}
\title{Return the dimension of a Codelink object.}
\description{
	Takes a Codelink object and returns the dimension (genes x samples).
}
\usage{
  \method{dim}{Codelink}(x)
}
\arguments{
  \item{x}{an object of class "Codelink".}
}
\value{A numeric vector with the dimensions.}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
\dontrun{
	data(codelink.example)
	dim(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
