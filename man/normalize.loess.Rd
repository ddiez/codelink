\name{normalize.loess}
\alias{normalize.loess}
\title{Normalized chips using cyclic loess.}
\description{
	Takes a matrix and apply cyclic loess normalization. It is based in
	normalize.loess from package affy but supports NA.
}
\usage{
	normalize.loess(mat, subset = sample(1:(dim(mat)[1]),
            min(c(5000,	nrow(mat)))), epsilon = 10^-2, maxit = 1,
            log.it = TRUE, verbose = TRUE, span = 2/3,
            family.loess = "symmetric")
}
\arguments{
  \item{mat}{a matrix with columns containing the values of the chips
    to normalize.}
  \item{subset}{a subset of the data to fit a loess to.}
  \item{epsilon}{a tolerance value (supposed to be a small value - used
    as a stopping criterium).}
  \item{maxit}{maximum number of iterations.}
  \item{log.it}{logical. If \code{TRUE} it takes the log2 of \code{mat}}
  \item{verbose}{logical. If \code{TRUE} displays current pair of chip being
    worked on.}
  \item{span}{parameter to be passed the function \code{\link[stats]{loess}}}
  \item{family.loess}{parameter to be passed the function
    \code{\link[stats]{loess}}. \code{"gaussian"} or \code{"symmetric"}
    are acceptable values for this parameter.}
}
\value{A matrix of normalized values.}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
#	mat <- matrix(sample(500), 100, 5)
#	mat <- normalize.loess(mat)
}
\keyword{documentation}
\keyword{utilities}
