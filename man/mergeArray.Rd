\name{mergeArray}
\alias{mergeArray}
\title{Merge Codelink Bioarrays Data}
\description{
  Merge data in a Codelink Object corresponding to same samples. Need a vector
indicating the classes and an optional vector indicating the labels of the mer-
ged samples.
}
\usage{
  mergeArray(object, class, names=NULL, method="mean",
             log.it=FALSE, merge.snr=TRUE)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{class}{a numeric vector indicating the classes.}
  \item{names}{an optional character vector indicating labels for each class.}
  \item{method}{the method used to summarize. Currently only "mean" supported.}
  \item{log.it}{logical; a logical indicating if log2 values should be returned.}
  \item{merge.snr}{logical; a logical indicating if SNR values should be merged.}
}
\value{An object of class "Codelink".}
\details{
}
\note{
}
\author{Diego Diez}
\seealso{
}
\examples{
#	data(codelink.example)
#	codelink.example <- bkgdCorrect(codelink.example)
#	codelink.example <- normalize(codelink.example, log.it=FALSE)
#	codelink.example <- mergeArray(codelink.example, class=c(1,1),
#					names="SAMPLE", log.it=TRUE)
}
\keyword{documentation}
\keyword{utilities}
