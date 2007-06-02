\name{normalize}
\alias{normalize}
\title{Normalization wrapper for Codelink objects.}
\description{
	Takes a Codelink object and applies normalization to intensity values.
}
\usage{
    normalize(object, method = "quantiles", log.it = TRUE, preserve = FALSE,
		weights = NULL, verbose = FALSE)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{method}{method to use in normalization.}
  \item{log.it}{logical; if data should be log2.}
  \item{preserve}{logical; if Ri slot should be preserved.}
  \item{weights}{weights vector for method CyclicLoess.}
  \item{verbose}{should informative output be printed.}
}
\value{A Codelink object with normalized intensity values.}
\details{
  Currently supported methods include "loess", "quantiles" and "median".  
  Median normalization is analogous to the default method applied for the 
  manufacturer in the Codelink software. Loess is a modified version of 
  CyclicLoess implemented in the affy package, allowing missing values and 
  weights. Quantile normalization uses the normalizeQuantiles() function in the 
  limma package.
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
#	data(codelink.example)
	# Background correction.
#	codelink.example <- bkgdCorrect(codelink.example,
#        method="half")
	# Normalization.
#	codelink.example <- normalize(codelink.example,
#        method="quantile")
}
\keyword{documentation}
\keyword{utilities}
