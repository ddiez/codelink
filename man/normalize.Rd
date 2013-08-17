\name{normalize}
\alias{normalize}
\alias{normalize,Codelink-method}
\alias{normalize,CodelinkSet-method}
\alias{codNormalize}
\alias{codNormalize,CodelinkSet-method}

\title{Normalization wrapper for Codelink/CodelinkSet objects.}
\description{
	Takes a CodelinkSet or Codelink object and applies normalization to intensity values.
}
\usage{

\S4method{normalize}{Codelink}(object, method="quantiles", log.it=TRUE, preserve=FALSE, weights=NULL, loess.method="fast")
	
\S4method{normalize}{CodelinkSet}(object, method="quantile", log.it=TRUE, weights=NULL, loess.method="fast")
    
\S4method{codNormalize}{CodelinkSet}(object, method="quantile", log.it=TRUE, weights=NULL, loess.method="fast")
}
\arguments{
  \item{object}{an object of class CodelinkSet or Codelink.}
  \item{...}{further arguments passed to the normalization functions.}
  \item{method}{method to use in normalization.}
  \item{log.it}{logical; if data should be log2.}
  \item{preserve}{logical; if Ri slot should be preserved.}
  \item{weights}{weights for method CyclicLoess.}
  \item{loess.method}{loess method to be used, default to 'fast'.}
}
\value{A CodelinkSet or Codelink object with normalized intensity values.}
\details{
  Currently supported methods include "median", "quantile" and "loess". Normalization functions are borrowed from the limma package.
  
  If a matrix of weights is passed for CyclicLoess normalization, a vector is constructed taking the lowest weight for each probe.
}
\author{Diego Diez}
\examples{
	data(codelink.exprset)
	# Background correction.
	codelink.exprset <- codCorrect(codelink.exprset,
        method = "half")
	# Normalization.
	codelink.exprset <- normalize(codelink.exprset,
        method = "quantile")
}
\keyword{documentation}
\keyword{utilities}
