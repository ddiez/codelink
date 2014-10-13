\name{averageProbes}
\alias{averageProbes}
\alias{averageProbes,CodelinkSet-method}
\title{averageProbes}
\description{
	Takes a CodelinkSet object and computes the average (mean) and sd of duplicated
	probes.
}
\usage{
	averageProbes(object, parallel = FALSE)
}
\arguments{
  \item{object}{an object of class "CodelinkSet".}
  \item{parallel}{whether to use a parallel version (requires package multicore).}
}
\details{

	This function will compute the mean() and sd() on each duplicated probe (i.e.
	identical probe id as for probeNames(). CodelinkSet objects use the per-array
	index to enable unique identifiers needed for ExpressionSet derived objects.
	Although the method probeNames() provides access to Codelink probe ids, this is
	inconvenient when dealing with other packages that make use of featureNames()
	to obtain probe ids and feed them  to the corresponding annotation package.

	In such cases CodelinkSet objects are not compatible with methods working on
	ExpressionSet objects. To avoid this limitation it is possible to construct a
	CodelinkSet object containing unique identifiers (CodelinkSetUnique class) by
	averaging the intensities of all replicated probes. This is done by computing
	the mean(). The sd() is also computed and stored in the slot sd.

	The current implementation takes a lot of time so a parallelized version of 
	lapply() may be used through the package 'multicore'. This is controlled by the
	argument 'parallel' which is FALSE by default.

	According to the authors it is not actually possible to use 'multicore' in a 
	GUI environment like the R.app Cocoa application in R, because this interferes 
	with the events loop. Therefore it is adviced to use option 'parallel=TRUE' in
	an R session running in a shell.

}
\author{Diego Diez}
\examples{
\dontrun{
    data(codelink.example)
	foo <- averageProbes(codelink.example)
}
}
\keyword{documentation}
\keyword{utilities}
