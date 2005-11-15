\name{normalize.Codelink}
\alias{normalize.Codelink}
\title{Normalization wrapper for Codelink objects.}
\description{
	Takes a Codelink object and applies normalization to intensity values.
}
\usage{
normalize.Codelink(object, method="quantiles", log.it=TRUE, type="ALL")
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{method}{method to use in normalization.}
	\item{log.it}{logical; if data should be log2.}
	\item{type}{which spots should be used.}
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
	codelink.example <- bkgdCorrect.Codelink(codelink.example,
        method="half")
	# Normalization.
	codelink.example <- normalize.Codelink(codelink.example,
        method="quantile")
}
\keyword{documentation}
\keyword{utilities}
