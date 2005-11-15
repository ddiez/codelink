\name{mergeArray.Codelink}
\alias{mergeArray.Codelink}
\title{Merge Codelink Bioarrays Data}
\description{
  Merge data in a Codelink Object corresponding to same samples. Need a vector
indicating the classes and an optional vector indicating the labels of the mer-
ged samples.
}
\usage{
	mergeArray.Codelink(object, group, names=NULL, method="mean", log.it=FALSE)
}
\arguments{
  \item{object}{an object of class "Codelink".}
  \item{group}{a numeric vector indicating the classes.}
  \item{names}{an optional character vector indicating labels for each class.}
  \item{method}{the method used to summarize. Currently only "mean" supported.}
  \item{log.it}{logical; a logical indicating if log2 values should be returned.}
}
\value{An object of class "Codelink".}
\details{
	TODO
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
	data(codelink.example)
	codelink.example <- bkgdCorrect.Codelink(codelink.example)
	codelink.example <- normalize.Codelink(codelink.example, log.it=FALSE)
	codelink.example <- mergeArray.Codelink(codelink.example,
	group=c(1,1), names="SAMPLE", log.it=TRUE)
}
\keyword{documentation}
\keyword{utilities}
