\name{normalize}
\alias{normalize}
\alias{normalize,Codelink-method}
\alias{normalize,CodelinkSet-method}
\alias{codNormalize}
\alias{codNormalize,CodelinkSet-method}

\title{Normalization wrapper for Codelink/CodelinkSet objects.}
\description{
	Takes a Codelink object and applies normalization to intensity values.
}
\usage{
    normalize(object, ...)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{...}{further arguments passed to the normalization functions.}
%  \item{method}{method to use in normalization.}
%  \item{log.it}{logical; if data should be log2.}
%  \item{preserve}{logical; if Ri slot should be preserved.}
%  \item{weights}{weights vector for method CyclicLoess.}
%  \item{verbose}{should informative output be printed.}
}
\value{A Codelink object with normalized intensity values.}
\details{
  Currently supported methods include "median", "quantile" and "loess". Normalization functions are borrowed from the limma package.
}
\author{Diego Diez}
\examples{
\dontrun{
	data(codelink.example)
	# Background correction.
	codelink.example <- bkgdCorrect(codelink.example,
        method = "half")
	# Normalization.
	codelink.example <- normalize(codelink.example,
        method = "quantile")
}
}
\keyword{documentation}
\keyword{utilities}
