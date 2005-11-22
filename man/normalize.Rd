\name{normalize}
\alias{normalize}
\title{Normalization wrapper for Codelink objects.}
\description{
	Takes a Codelink object and applies normalization to intensity values.
}
\usage{
    normalize(object, method="quantiles", log.it=TRUE)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{method}{method to use in normalization.}
	\item{log.it}{logical; if data should be log2.}
}
\value{A Codelink object with normalized intensity values.}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
	data(codelink.example)
	# Background correction.
	codelink.example <- bkgdCorrect(codelink.example,
        method="half")
	# Normalization.
	codelink.example <- normalize(codelink.example,
        method="quantile")
}
\keyword{documentation}
\keyword{utilities}
