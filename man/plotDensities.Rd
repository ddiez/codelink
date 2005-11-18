\name{plotDensities}
\alias{plotDensities}
\title{Plot Densities}
\description{
	Takes a Codelink object and plot the distributions of intensities.
}
\usage{
	plotDensities(object, subset=c(1:dim(object)[2]), title=NULL, legend.cex=1)
}
\arguments{
	\item{object}{an object of class "Codelink".}
	\item{subset}{subset of arrays to be plotted (default: all).}
	\item{title}{title of the plot.}
	\item{legend.cex}{font factor use in legend to fit.}
}
\value{
}
\details{
}
\note{
}
\author{Diego Díez Ruiz}
\seealso{
}
\examples{
    data(codelink.example)
    plotDensities(codelink.example)
}
\keyword{documentation}
\keyword{utilities}
